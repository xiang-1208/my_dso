#pragma once
#include "util/NumType.h"
#include "util/globalutil.h"
#include "FullSystem/PixelSelector2.h"
#include "util/nanoflann.h"
#include "IOWrapper/Output3DWrapper.h"
#include "OptimizationBackend/MatrixAccumulators.h"
#include "util/globalFuncs.h"

using namespace dso;

struct Pnt
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

	// index in jacobian. never changes (actually, there is no reason why).
	float u,v;

	// idepth / isgood / energy during optimization.
	float idepth;				//!< 该点对应参考帧的逆深度
	float idepth_new;			//!< 该点在新的一帧(当前帧)上的逆深度
	float iR;					//!< 逆深度的期望值

	bool isGood;				//!< 点在新图像内, 相机前, 像素值有穷则好
	bool isGood_new;

	Vec2f energy;				//!< [0]残差的平方, [1]正则化项(逆深度减一的平方)
	Vec2f energy_new;			//!< 迭代计算的新的能量

	float lastHessian;			//!< 逆深度的Hessian, 即协方差, dd*dd
	float lastHessian_new;		//!< 新一次迭代的协方差

	// max stepsize for idepth (corresponding to max. movement in pixel-space).
	float maxstep;				//!< 逆深度增加的最大步长

	float my_type; 				//!< 第0层提取是1, 2, 4, 对应d, 2d, 4d, 其它层是1
	float outlierTH; 			//!< 外点阈值

	// idx (x+y*w) of up to 10 nearest points in pixel space.
	int neighbours[10];			//!< 图像中离该点最近的10个点
	float neighboursDist[10];   //!< 最近10个点的距离

	// idx (x+y*w) of closest point one pyramid level above.
	int parent;		  			//!< 上一层中该点的父节点 (距离最近的)的id
	float parentDist;			//!< 上一层中与父节点的距离
};

class CoarseInitializer
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	CoarseInitializer(int ww, int hh);  
    void setFirst(CalibHessian* HCalib, FrameHessian* newFrameHessian);
	bool trackFrame(FrameHessian* newFrameHessian, std::vector<IOWrap::Output3DWrapper*> &wraps);


	Pnt* points[PYR_LEVELS];
	int numPoints[PYR_LEVELS];
	AffLight thisToNext_aff;
	SE3 thisToNext;

    int frameID;

protected:
    void makeK(CalibHessian* HCalib);
	void makeNN();
	void propagateDown(int srcLvl);
	void optReg(int lvl);
	void resetPoints(int lvl);
	Vec3f calcResAndGS(
		int lvl, Mat88f &H_out, Vec8f &b_out,
		Mat88f &H_out_sc, Vec8f &b_out_sc,
		const SE3 &refToNew, AffLight refToNew_aff,
		bool plot);

	Mat33 K[PYR_LEVELS];
	Mat33 Ki[PYR_LEVELS];
	double fx[PYR_LEVELS];
	double fy[PYR_LEVELS];
	double fxi[PYR_LEVELS];
	double fyi[PYR_LEVELS];
	double cx[PYR_LEVELS];
	double cy[PYR_LEVELS];
	double cxi[PYR_LEVELS];
	double cyi[PYR_LEVELS];
	int w[PYR_LEVELS];
	int h[PYR_LEVELS];

	//? 这几个参数很迷
	float alphaK;					//!< 2.5*2.5
	float alphaW;					//!< 150*150
	float regWeight;				//!< 对逆深度的加权值, 0.8
	float couplingWeight;			//!< 1

    FrameHessian* firstFrame;
	FrameHessian* newFrame;			//!< track中新加入的帧
    
	bool snapped;					//!< 是否尺度收敛 (暂定)
	int snappedAt;					//!< 尺度收敛在第几帧

	Vec3f dGrads[PYR_LEVELS];

	//* 9维向量, 乘积获得9*9矩阵, 并做的累加器
	Accumulator9 acc9;			//!< Hessian 矩阵

	// temporary buffers for H and b.
	Vec10f* JbBuffer;			//!< 用来计算Schur的 0-7: sum(dd * dp). 8: sum(res*dd). 9: 1/(1+sum(dd*dd))=inverse hessian entry.
	Vec10f* JbBuffer_new;		//!< 跌待更新后新的值
};

//* 作为 KDTreeSingleIndexAdaptor 类的第二个模板参数必须给出, 包括下面的接口
struct FLANNPointcloud
{
    inline FLANNPointcloud() {num=0; points=0;}
    inline FLANNPointcloud(int n, Pnt* p) :  num(n), points(p) {}
	int num;
	Pnt* points;
	// 返回数据点的数目
	inline size_t kdtree_get_point_count() const { return num; }
	// 使用L2度量时使用, 返回向量p1, 到第idx_p2个数据点的欧氏距离
	inline float kdtree_distance(const float *p1, const size_t idx_p2,size_t /*size*/) const
	{
		const float d0=p1[0]-points[idx_p2].u;
		const float d1=p1[1]-points[idx_p2].v;
		return d0*d0+d1*d1;
	}
	// 返回第idx个点的第dim维数据
	inline float kdtree_get_pt(const size_t idx, int dim) const
	{
		if (dim==0) return points[idx].u;
		else return points[idx].v;
	}
	// 可选计算bounding box
	// false 表示默认
	// true 本函数应该返回bb
	template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const { return false; }
};