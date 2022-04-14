#pragma once
#include "Eigen/Core"
#include "sophus/sim3.hpp"
#include "sophus/se3.hpp"

#define CPARS 4
#define MAX_RES_PER_POINT 8

typedef Sophus::SE3d SE3;

#define MAX_RES_PER_POINT 8

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> MatXX;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
typedef Eigen::Matrix<float, MAX_RES_PER_POINT,1> VecNRf;
typedef Eigen::Matrix<double,3,3> Mat33;
typedef Eigen::Matrix<double,4,2> Mat42;
typedef Eigen::Matrix<double,6,6> Mat66;
typedef Eigen::Matrix<double,8,8> Mat88;
typedef Eigen::Matrix<float,2,2> Mat22f;
typedef Eigen::Matrix<float,1,8> Mat18f;
typedef Eigen::Matrix<float,8,8> Mat88f;
typedef Eigen::Matrix<float,3,3> Mat33f;
typedef Eigen::Matrix<float,9,9> Mat99f;

typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;

typedef Eigen::Matrix<double,2,1> Vec2;
typedef Eigen::Matrix<float,9,1> Vec9f;
typedef Eigen::Matrix<double,CPARS,1> VecC;
typedef Eigen::Matrix<float,CPARS,1> VecCf;
typedef Eigen::Matrix<float,3,1> Vec3f;
typedef Eigen::Matrix<float,2,1> Vec2f;
typedef Eigen::Matrix<float,8,1> Vec8f;
typedef Eigen::Matrix<float,10,1> Vec10f;
typedef Eigen::Matrix<double,3,1> Vec3;
typedef Eigen::Matrix<double,6,1> Vec6;
typedef Eigen::Matrix<double,8,1> Vec8;
typedef Eigen::Matrix<double,10,1> Vec10;
typedef Eigen::Matrix<unsigned char,3,1> Vec3b;

struct AffLight
{
	AffLight(double a_, double b_) : a(a_), b(b_) {};
	AffLight() : a(0), b(0) {};    

    double a,b; // I_frame = exp(a)*I_global + b. // I_global = exp(-a)*(I_frame - b).

	Vec2 vec()
	{
		return Vec2(a,b);
	}

	/********************************
	 * @ function: 把光度仿射变换转化为, 能量函数中的整体光度仿射系数
	 * @
	 * @ param: 	exposureF		参考帧曝光时间
	 * @ 			exposureT		目标帧曝光时间
	 * @			g2F				参考帧光度仿射系数
	 * @			g2T				目标帧光度仿射系数
	 * @ note:	注意这里面的a,b之间的差别
	 *******************************/	
	static Vec2 fromToVecExposure(float exposureF, float exposureT, AffLight g2F, AffLight g2T)
	{
		// 没有曝光时间标定, 置1
		if(exposureF==0 || exposureT==0)
		{
			exposureT = exposureF = 1;
			//printf("got exposure value of 0! please choose the correct model.\n");
			//assert(setting_brightnessTransferFunc < 2);
		}

		// 论文公式(4), 光度系数放一起
		double a = exp(g2T.a-g2F.a) * exposureT / exposureF;
		double b = g2T.b - a*g2F.b;
		return Vec2(a,b);
	}
};