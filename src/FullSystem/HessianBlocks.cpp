#include "HessianBlocks.h"
#include "OptimizationBackend/EnergyFunctionalStructs.h"

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