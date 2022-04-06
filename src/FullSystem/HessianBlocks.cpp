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