#pragma once
#include "Eigen/Core"
#include "sophus/sim3.hpp"
#include "sophus/se3.hpp"

#define CPARS 4

typedef Sophus::SE3d SE3;

#define MAX_RES_PER_POINT 8

typedef Eigen::Matrix<double,Eigen::Dynamic,1> VecX;
typedef Eigen::Matrix<double,3,3> Mat33;
typedef Eigen::Matrix<float,8,8> Mat88f;

typedef Eigen::Matrix<double,CPARS,1> VecC;
typedef Eigen::Matrix<float,CPARS,1> VecCf;
typedef Eigen::Matrix<float,3,1> Vec3f;
typedef Eigen::Matrix<float,2,1> Vec2f;
typedef Eigen::Matrix<float,8,1> Vec8f;
typedef Eigen::Matrix<unsigned char,3,1> Vec3b;

struct AffLight
{
	AffLight(double a_, double b_) : a(a_), b(b_) {};
	AffLight() : a(0), b(0) {};    

    double a,b; // I_frame = exp(a)*I_global + b. // I_global = exp(-a)*(I_frame - b).
};