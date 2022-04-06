#include "ImmaturePoint.h"
#include "HessianBlocks.h"

ImmaturePoint::ImmaturePoint(int u_, int v_, FrameHessian* host_, float type, CalibHessian* HCalib)
: u(u_), v(v_), host(host_), my_type(type), idepth_min(0), idepth_max(NAN), lastTraceStatus(IPS_UNINITIALIZED)
{
    gradH.setZero();

    for(int idx=0;idx<patternNum;idx++)
    {
		int dx = patternP[idx][0];
		int dy = patternP[idx][1];    

		// 由于+0.5导致积分, 插值得到值3个 [像素值, dx, dy]
        Vec3f ptc = getInterpolatedElement33BiLin(host->dI, u+dx, v+dy,wG[0]);      

        color[idx] = ptc[0];  
        if(!std::isfinite(color[idx])) {energyTH=NAN; return;}

		// 梯度矩阵[dx*2, dxdy; dydx, dy^2]
		gradH += ptc.tail<2>()  * ptc.tail<2>().transpose();    
		//! 点的权重 c^2 / ( c^2 + ||grad||^2 )
		weights[idx] = sqrtf(setting_outlierTHSumComponent / (setting_outlierTHSumComponent + ptc.tail<2>().squaredNorm()));         
    }
	energyTH = patternNum*setting_outlierTH;
	energyTH *= setting_overallEnergyTHWeight*setting_overallEnergyTHWeight;

	idepth_GT=0;
	quality=10000;
}