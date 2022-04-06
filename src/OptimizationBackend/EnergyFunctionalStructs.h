#pragma once

#include "util/NumType.h"
#include "vector"

class FrameHessian; //前向声明
class EFResidual;
class EFPoint;
class PointHessian;

class EFFrame
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;   
    EFFrame(FrameHessian* d) : data(d)
	{
		takeData();  
	}

    void takeData();

	//! 位姿 0-5, 光度ab 6-7
	Vec8 prior;					//!< 位姿只有第一帧有先验 prior hessian (diagonal)
	Vec8 delta_prior;			//!< 相对于先验的增量   // = state-state_prior (E_prior = (delta_prior)' * diag(prior) * (delta_prior)
	Vec8 delta;					//!< 相对于线性化点位姿, 光度的增量  // state - state_zero.

    FrameHessian* data;				//!< 对应FrameHessian数据 
	std::vector<EFPoint*> points;	//!< 帧上所有点

	int idx;	//!< 在能量函数中帧id // idx in frames.
	int frameID;	//!< 所有历史帧ID 

};

class EFPoint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	// contains all residuals.
	std::vector<EFResidual*> residualsAll;	//!< 该点的所有残差

	float deltaF;		//!< 当前逆深度和线性化处的差, 没有使用FEJ, 就是0

	PointHessian* data;	//!< PointHessian数据


};

class EFResidual
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EFPoint* point;					//!< 残差点

	EFFrame* host;					//!< 主
	EFFrame* target;				//!< 目标
	int hostIDX, targetIDX;  		//!< 残差对应的 host 和 Target ID号
};