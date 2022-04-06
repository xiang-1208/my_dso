#pragma once
#include <Eigen/Core>
#include "util/globalFuncs.h"
#include "util/NumType.h"

class FrameHessian;
class CalibHessian;


enum ImmaturePointStatus {
	IPS_GOOD=0,					// traced well and good
	// 搜索区间超出图像, 尺度变化太大, 两次残差都大于阈值, 不再搜索
	IPS_OOB,					// OOB: end tracking & marginalize!
	// 第一次残差大于阈值, 外点
	IPS_OUTLIER,				// energy too high: if happens again: outlier!
	// 搜索区间太短，但是没激活
	IPS_SKIPPED,				// traced well and good (but not actually traced).
	// 梯度和极线夹角太大
	IPS_BADCONDITION,			// not traced because of bad condition.
	IPS_UNINITIALIZED};			// not even traced once.



class ImmaturePoint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;  
    float color[MAX_RES_PER_POINT];		//!< 原图上pattern上对应的像素值
    float weights[MAX_RES_PER_POINT];	//!< 原图上pattern对应的权重(与梯度成反比)

    Mat22f gradH;				//!< 图像梯度hessian矩阵

    float u,v;					//!< host里的像素坐标  
    FrameHessian* host;
	float my_type;
    float energyTH;

	float idepth_min;			//!< 逆深度范围
	float idepth_max;

    float idepth_GT;
    float quality;				//!< 第二误差/第一误差 作为搜索质量, 越大越好

    ImmaturePointStatus lastTraceStatus;		//!< 上一次跟踪状态

    ImmaturePoint(int u_, int v_, FrameHessian* host_, float type, CalibHessian* HCalib);
};