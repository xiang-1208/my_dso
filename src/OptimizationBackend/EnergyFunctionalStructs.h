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

enum EFPointStatus {PS_GOOD=0, PS_MARGINALIZE, PS_DROP};

class EFPoint
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	EFPoint(PointHessian* d, EFFrame* host_) : data(d),host(host_)
	{
		takeData();
		stateFlag=EFPointStatus::PS_GOOD;
	}

	void takeData();

	// contains all residuals.
	std::vector<EFResidual*> residualsAll;	//!< 该点的所有残差

	float priorF;		//!< 逆深度先验信息矩阵, 初始化之后的有
	float deltaF;		//!< 当前逆深度和线性化处的差, 没有使用FEJ, 就是0

	EFFrame* host;

	int idxInPoints;	//!< 当前点在EFFrame中id

	PointHessian* data;	//!< PointHessian数据

	EFPointStatus stateFlag; //!< 点的状态
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