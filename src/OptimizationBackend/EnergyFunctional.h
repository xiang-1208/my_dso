#pragma once


#include "util/NumType.h"
#include "FullSystem/HessianBlocks.h"
#include "EnergyFunctionalStructs.h"
#include "map"

extern bool EFAdjointsValid;		//!< 是否设置状态伴随矩阵
extern bool EFIndicesValid;		//!< 是否设置frame, point, res的ID
extern bool EFDeltaValid;			//!< 是否设置状态增量值

class EnergyFunctional
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EnergyFunctional();
    
    void setAdjointsF(CalibHessian* Hcalib);
	void setDeltaF(CalibHessian* HCalib);

    EFFrame* insertFrame(FrameHessian* fh, CalibHessian* Hcalib);
    void makeIDX();


    std::vector<EFFrame*> frames;  		//!< 能量函数中的帧
    int nPoints, nFrames, nResiduals;	//!< EFPoint的数目, EFframe关键帧数, 残差数

	MatXX HM;					//!< 优化的Hessian矩阵, 边缘化掉逆深度
	VecX bM;					//!< 优化的Jr项, 边缘化掉逆深度

	EFPoint* insertPoint(PointHessian* ph);

	std::map<uint64_t, // 历史ID
	    Eigen::Vector2i,
	    std::less<uint64_t>, 
	    Eigen::aligned_allocator<std::pair<const uint64_t, Eigen::Vector2i>> // 64位对齐
	    > connectivityMap; 			//!< 关键帧之间的连接关系, first: 前32表示host ID, 后32位表示target ID; second:数目 [0] 普通的, [1] 边缘化的
protected: 
	Mat18f* adHTdeltaF;				//!< host和target之间位姿的增量, 一共帧数×帧数个

	Mat88* adHost; 					//!< 伴随矩阵, double
	Mat88* adTarget;  

	Mat88f* adHostF;				//!< 伴随矩阵, float
	Mat88f* adTargetF;   

	VecCf cDeltaF;									//!< 相机内参增量
    VecC cPrior;		//!< setting_initialCalibHessian 信息矩阵 
    VecCf cPriorF;

    std::vector<EFPoint*> allPoints;
};