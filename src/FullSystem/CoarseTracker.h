#pragma once

 
#include "util/NumType.h"
#include "vector"
#include <math.h>
#include "util/globalutil.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "IOWrapper/Output3DWrapper.h"
#include "FullSystem/HessianBlocks.h"
#include <algorithm>
#include <stdlib.h>

class CoarseTracker
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    CoarseTracker(int w, int h);

	FrameHessian* lastRef;			//!< 参考帧
	AffLight lastRef_aff_g2l;
	FrameHessian* newFrame;			//!< 新来的一帧
	int refFrameID;					//!< 参考帧id

	bool debugPrint, debugPlot;

	Mat33f K[PYR_LEVELS];
	Mat33f Ki[PYR_LEVELS];
	float fx[PYR_LEVELS];
	float fy[PYR_LEVELS];
	float fxi[PYR_LEVELS];
	float fyi[PYR_LEVELS];
	float cx[PYR_LEVELS];
	float cy[PYR_LEVELS];
	float cxi[PYR_LEVELS];
	float cyi[PYR_LEVELS];
	int w[PYR_LEVELS];
	int h[PYR_LEVELS];

private:
	float* idepth[PYR_LEVELS];
	float* weightSums[PYR_LEVELS];
	float* weightSums_bak[PYR_LEVELS];

	// pc buffers
	float* pc_u[PYR_LEVELS];				//!< 每层上的有逆深度点的坐标x
	float* pc_v[PYR_LEVELS];				//!< 每层上的有逆深度点的坐标y
	float* pc_idepth[PYR_LEVELS];			//!< 每层上点的逆深度
	float* pc_color[PYR_LEVELS];			//!< 每层上点的颜色值
	int pc_n[PYR_LEVELS];					//!< 每层上点的个数

	// warped buffers
	float* buf_warped_idepth;				//!< 投影得到的点的逆深度
	float* buf_warped_u;					//!< 投影得到的归一化坐标
	float* buf_warped_v;					//!< 同上
	float* buf_warped_dx;					//!< 投影点的图像梯度
	float* buf_warped_dy;					//!< 同上
	float* buf_warped_residual;				//!< 投影得到的残差
	float* buf_warped_weight;				//!< 投影的huber函数权重
	float* buf_warped_refColor;				//!< 投影点参考帧上的灰度值
	int buf_warped_n;						//!< 投影点的个数  

	std::vector<float*> ptrToDelete;				//!< 所有的申请的内存指针, 用于析构删除  
};