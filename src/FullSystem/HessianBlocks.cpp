#include "HessianBlocks.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"

//@ 计算优化前和优化后的相对位姿, 相对光度变化, 及中间变量
void FrameFramePrecalc::set(FrameHessian* host, FrameHessian* target, CalibHessian* HCalib )
{
	this->host = host;    // 这个是赋值, 计数会增加, 不是拷贝
	this->target = target;

	//? 实在不懂leftToleft_0这个名字怎么个含义
	// 优化前host target间位姿变换
	SE3 leftToLeft_0 = target->get_worldToCam_evalPT() * host->get_worldToCam_evalPT().inverse();
	PRE_RTll_0 = (leftToLeft_0.rotationMatrix()).cast<float>();
	PRE_tTll_0 = (leftToLeft_0.translation()).cast<float>();   

	// 优化后host到target间位姿变换
	SE3 leftToLeft = target->PRE_worldToCam * host->PRE_camToWorld;
	PRE_RTll = (leftToLeft.rotationMatrix()).cast<float>();
	PRE_tTll = (leftToLeft.translation()).cast<float>();
	distanceLL = leftToLeft.translation().norm();   

	// 乘上内参, 中间量?
	Mat33f K = Mat33f::Zero();
	K(0,0) = HCalib->fxl();
	K(1,1) = HCalib->fyl();
	K(0,2) = HCalib->cxl();
	K(1,2) = HCalib->cyl();
	K(2,2) = 1;
	PRE_KRKiTll = K * PRE_RTll * K.inverse();
	PRE_RKiTll = PRE_RTll * K.inverse();
	PRE_KtTll = K * PRE_tTll;

	// 光度仿射值
	PRE_aff_mode = AffLight::fromToVecExposure(host->ab_exposure, target->ab_exposure, host->aff_g2l(), target->aff_g2l()).cast<float>();
	PRE_b0_mode = host->aff_g2l_0().b;  
}

//@ 设置固定线性化点位置的状态
//TODO 后面求nullspaces地方没看懂, 回头再看<2019.09.18> 数学原理是啥?
void FrameHessian::setStateZero(const Vec10 &state_zero)
{
	//! 前六维位姿必须是0
	assert(state_zero.head<6>().squaredNorm() < 1e-20);

	this->state_zero = state_zero;

	//! 感觉这个nullspaces_pose就是 Adj_T
	//! Exp(Adj_T*zeta)=T*Exp(zeta)*T^{-1}
	// 全局转为局部的，左乘边右乘
	//! T_c_w * delta_T_g * T_c_w_inv = delta_T_l
	//TODO 这个是数值求导的方法么???
	for(int i=0;i<6;i++)
	{
		Vec6 eps; eps.setZero(); eps[i] = 1e-3;
		SE3 EepsP = Sophus::SE3::exp(eps);
		SE3 EepsM = Sophus::SE3::exp(-eps);
		SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT() * EepsP) * get_worldToCam_evalPT().inverse();
		SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT() * EepsM) * get_worldToCam_evalPT().inverse();
		nullspaces_pose.col(i) = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);
	}	
	//nullspaces_pose.topRows<3>() *= SCALE_XI_TRANS_INVERSE;
	//nullspaces_pose.bottomRows<3>() *= SCALE_XI_ROT_INVERSE;

	//? rethink
	// scale change
	SE3 w2c_leftEps_P_x0 = (get_worldToCam_evalPT());
	w2c_leftEps_P_x0.translation() *= 1.00001;
	w2c_leftEps_P_x0 = w2c_leftEps_P_x0 * get_worldToCam_evalPT().inverse();
	SE3 w2c_leftEps_M_x0 = (get_worldToCam_evalPT());
	w2c_leftEps_M_x0.translation() /= 1.00001;
	w2c_leftEps_M_x0 = w2c_leftEps_M_x0 * get_worldToCam_evalPT().inverse();
	nullspaces_scale = (w2c_leftEps_P_x0.log() - w2c_leftEps_M_x0.log())/(2e-3);


	nullspaces_affine.setZero();
	nullspaces_affine.topLeftCorner<2,1>()  = Vec2(1,0);
	assert(ab_exposure > 0);
	nullspaces_affine.topRightCorner<2,1>() = Vec2(0, expf(aff_g2l_0().a)*ab_exposure);
}

void FrameHessian::makeImages(float* color, CalibHessian* HCalib)
{
    for(int i=0;i<pyrLevelsUsed;i++)
    {
        dIp[i] = new Eigen::Vector3f[wG[i]*hG[i]];
        absSquaredGrad[i] = new float[wG[i]*hG[i]];
    }
    dI = dIp[0];

    int w=wG[0];
    int h=hG[0];
	for(int i=0;i<w*h;i++)
    {
		dI[i][0] = color[i];
        // std::cout << dI[i][0] << "    ";
    }

    for(int lvl=0; lvl<pyrLevelsUsed; lvl++)
    {
        int wl = wG[lvl], hl = hG[lvl];
        Eigen::Vector3f* dI_l = dIp[lvl];

        float* dabs_l = absSquaredGrad[lvl];
        if(lvl > 0)
        {
            int lvlm1 = lvl-1;
            int wlm1 = wG[lvlm1];
            Eigen::Vector3f* dI_lm = dIp[lvlm1];

			for(int y=0;y<hl;y++)
				for(int x=0;x<wl;x++)
                {
                    dI_l[x+y*wl][0] = 0.25f * (dI_lm[2*x   + 2*y*wlm1][0] +
												dI_lm[2*x+1 + 2*y*wlm1][0] +
												dI_lm[2*x   + 2*y*wlm1+wlm1][0] +
												dI_lm[2*x+1 + 2*y*wlm1+wlm1][0]);
                }            
        }

        for(int idx=wl;idx < wl*(hl-1);idx++)
        {
            float dx = 0.5f*(dI_l[idx+1][0] - dI_l[idx-1][0]);
            float dy = 0.5f*(dI_l[idx+wl][0] - dI_l[idx-wl][0]);

			if(!std::isfinite(dx)) dx=0;
			if(!std::isfinite(dy)) dy=0;

			dI_l[idx][1] = dx;
			dI_l[idx][2] = dy;      

            dabs_l[idx] = dx*dx+dy*dy;    

			if(setting_gammaWeightsPixelSelect==1 && HCalib!=0)
			{
				float gw = HCalib->getBGradOnly((float)(dI_l[idx][0]));
				dabs_l[idx] *= gw*gw;	// convert to gradient of original color space (before removing response).
			}              
        }
    }
}

PointHessian::PointHessian(const ImmaturePoint* const rawPoint, CalibHessian* Hcalib)
{
	energyTH = rawPoint->energyTH;
	host = rawPoint->host; // 主帧
	hasDepthPrior=false;

	// set static values & initialization.
	u = rawPoint->u;
	v = rawPoint->v;	
	assert(std::isfinite(rawPoint->idepth_max));
	//idepth_init = rawPoint->idepth_GT;

	int n = patternNum;
	memcpy(color, rawPoint->color, sizeof(float)*n);// 一个点对应8个像素
	memcpy(weights, rawPoint->weights, sizeof(float)*n);

	efPoint=0; // 指针=0
}

