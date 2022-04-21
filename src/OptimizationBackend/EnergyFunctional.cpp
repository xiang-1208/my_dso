#include "EnergyFunctional.h"

bool EFAdjointsValid = false;		//!< 是否设置状态伴随矩阵
bool EFIndicesValid = false;		//!< 是否设置frame, point, res的ID
bool EFDeltaValid = false;			//!< 是否设置状态增量值

EnergyFunctional::EnergyFunctional()
{
	adHost=0;
	adTarget=0;

	adHostF=0;
	adTargetF=0;
	adHTdeltaF=0;

	HM = MatXX::Zero(CPARS,CPARS); // 初始的, 后面增加frame改变
	bM = VecX::Zero(CPARS);

    nFrames = nResiduals = nPoints = 0;
}


//@ 向能量函数中增加一帧, 进行的操作: 改变正规方程, 重新排ID, 共视关系
EFFrame* EnergyFunctional::insertFrame(FrameHessian* fh, CalibHessian* Hcalib)
{
	// 建立优化用的能量函数帧. 并加进能量函数frames中
	EFFrame* eff = new EFFrame(fh);    
    eff->idx = frames.size(); 
    frames.push_back(eff);

    nFrames++;
    fh->efFrame = eff; // FrameHessian 指向能量函数帧

    assert(HM.cols() == 8*nFrames+CPARS-8);  // 边缘化掉一帧, 缺8个
	// 一个帧8个参数 + 相机内参
	bM.conservativeResize(8*nFrames+CPARS);
	HM.conservativeResize(8*nFrames+CPARS,8*nFrames+CPARS);
	// 新帧的块为0
	bM.tail<8>().setZero();
	HM.rightCols<8>().setZero();
	HM.bottomRows<8>().setZero();

    EFIndicesValid = false;
	EFAdjointsValid=false;
	EFDeltaValid=false;

    setAdjointsF(Hcalib); 	// 设置伴随矩阵
	makeIDX();				// 设置ID

	for(EFFrame* fh2 : frames)
	{
		// 前32位是host帧的历史ID, 后32位是Target的历史ID
        connectivityMap[(((uint64_t)eff->frameID) << 32) + ((uint64_t)fh2->frameID)] = Eigen::Vector2i(0,0);
		if(fh2 != eff)
            connectivityMap[(((uint64_t)fh2->frameID) << 32) + ((uint64_t)eff->frameID)] = Eigen::Vector2i(0,0);
	}

	return eff;
}

//@ 计算各种状态的相对量的增量
void EnergyFunctional::setDeltaF(CalibHessian* HCalib)
{
	if(adHTdeltaF != 0) delete[] adHTdeltaF;
	adHTdeltaF = new Mat18f[nFrames*nFrames];
	std::cout << "!!!!!!!!!!!!!!!!!!!!!!!" <<std::endl;	

	for(int h=0;h<nFrames;h++)
		for(int t=0;t<nFrames;t++)
		{
			int idx = h+t*nFrames;
			//! delta_th = Adj * delta_t or delta_th = Adj * delta_h
			// 加一起应该是, 两帧之间位姿变换的增量, 因为h变一点, t变一点
			adHTdeltaF[idx] = frames[h]->data->get_state_minus_stateZero().head<8>().cast<float>().transpose() * adHostF[idx]
					        +frames[t]->data->get_state_minus_stateZero().head<8>().cast<float>().transpose() * adTargetF[idx];						
		}	

	cDeltaF = HCalib->value_minus_value_zero.cast<float>(); // 相机内参增量

	for(EFFrame* f : frames)
	{
		f->delta = f->data->get_state_minus_stateZero().head<8>();  // 帧位姿增量
		f->delta_prior = (f->data->get_state() - f->data->getPriorZero()).head<8>(); // 先验增量

		for(EFPoint* p : f->points)
			p->deltaF = p->data->idepth - p->data->idepth_zero; // 逆深度的增量
	}

	EFDeltaValid = true;	
}

//@ 向能量函数中插入一个点, 放入对应的EFframe
EFPoint* EnergyFunctional::insertPoint(PointHessian* ph)
{
	EFPoint* efp = 0;
	efp = new EFPoint(ph, ph->host->efFrame);
	efp->idxInPoints = ph->host->efFrame->points.size();
	ph->host->efFrame->points.push_back(efp);

	nPoints++;
	ph->efPoint = efp;

	EFIndicesValid = false; // 有插入需要重新梳理残差的ID

	return efp;
}

//@ 设置EFFrame, EFPoint, EFResidual对应的 ID 号
void EnergyFunctional::makeIDX()
{
	// 重新赋值ID
	for(unsigned int idx=0;idx<frames.size();idx++)
		frames[idx]->idx = idx;	

	allPoints.clear();

	for(EFFrame* f : frames)
		for(EFPoint* p : f->points) 
		{
			allPoints.push_back(p);
			// 残差的ID号
			for(EFResidual* r : p->residualsAll)
			{
				r->hostIDX = r->host->idx;  // EFFrame的idx
				r->targetIDX = r->target->idx;
			}
		}

	EFIndicesValid=true;
}

void EnergyFunctional::setAdjointsF(CalibHessian* Hcalib)
{
	if(adHost != 0) delete[] adHost;
	if(adTarget != 0) delete[] adTarget;
	adHost = new Mat88[nFrames*nFrames];
	adTarget = new Mat88[nFrames*nFrames];   

	for(int h=0;h<nFrames;h++) // 主帧
		for(int t=0;t<nFrames;t++) // 目标帧
        {
			FrameHessian* host = frames[h]->data;
			FrameHessian* target = frames[t]->data;   

            SE3 hostToTarget = target->get_worldToCam_evalPT() * host->get_worldToCam_evalPT().inverse();     

			Mat88 AH = Mat88::Identity();
			Mat88 AT = Mat88::Identity();

			// 见笔记推导吧, 或者https://www.cnblogs.com/JingeTU/p/9077372.html
			AH.topLeftCorner<6,6>() = -hostToTarget.Adj().transpose();
			AT.topLeftCorner<6,6>() = Mat66::Identity();

			// 光度参数, 合并项对参数求导
			//! E = Ij - tj*exp(aj) / ti*exp(ai) * Ii - (bj - tj*exp(aj) / ti*exp(ai) * bi)
			//! a = - tj*exp(aj) / ti*exp(ai),  b = - (bj - tj*exp(aj) / ti*exp(ai) * bi) 
			Vec2f affLL = AffLight::fromToVecExposure(host->ab_exposure,target->ab_exposure,host->aff_g2l_0(),target->aff_g2l_0()).cast<float>();
			AT(6,6) = -affLL[0]; //! a'(aj)
			AH(6,6) = affLL[0];	 //! a'(ai)
			AT(7,7) = -1;		 //! b'(bj)
			AH(7,7) = affLL[0];	 //! b'(bi)

			AH.block<3,8>(0,0) *= SCALE_XI_TRANS;
			AH.block<3,8>(3,0) *= SCALE_XI_ROT;
			AH.block<1,8>(6,0) *= SCALE_A;
			AH.block<1,8>(7,0) *= SCALE_B;
			AT.block<3,8>(0,0) *= SCALE_XI_TRANS;
			AT.block<3,8>(3,0) *= SCALE_XI_ROT;
			AT.block<1,8>(6,0) *= SCALE_A;    //? 已经是乘过的, 怎么又乘一遍
			AT.block<1,8>(7,0) *= SCALE_B;

			adHost[h+t*nFrames] = AH;
			adTarget[h+t*nFrames] = AT;			
        } 
	cPrior = VecC::Constant(setting_initialCalibHessian); // 常数矩阵

	// float型
	if(adHostF != 0) delete[] adHostF;
	if(adTargetF != 0) delete[] adTargetF;
	adHostF = new Mat88f[nFrames*nFrames];
	adTargetF = new Mat88f[nFrames*nFrames];	

	for(int h=0;h<nFrames;h++)
		for(int t=0;t<nFrames;t++)
		{
			adHostF[h+t*nFrames] = adHost[h+t*nFrames].cast<float>();
			adTargetF[h+t*nFrames] = adTarget[h+t*nFrames].cast<float>();
		}

	cPriorF = cPrior.cast<float>();


	EFAdjointsValid = true;
}