#pragma once
#include <Eigen/Core>
#include "FullSystem/FrameShell.h"
#include "util/globalutil.h"

using namespace dso;

struct CalibHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static int instanceCounter;

    float Binv[256];
    float B[256];

    VecCf value_scaledf;

    inline ~CalibHessian() {instanceCounter--;}
	inline CalibHessian()
	{
		VecC initial_value = VecC::Zero();
		initial_value[0] = fxG[0];
		initial_value[1] = fyG[0];
		initial_value[2] = cxG[0];
		initial_value[3] = cyG[0];

		// setValueScaled(initial_value);
		// value_zero = value;
		// value_minus_value_zero.setZero();

		instanceCounter++;
		for(int i=0;i<256;i++)
			Binv[i] = B[i] = i;		// set gamma function to identity
	};

    EIGEN_STRONG_INLINE float getBGradOnly(float color)
	{
		int c = color+0.5f;
		if(c<5) c=5;
		if(c>250) c=250;
		return B[c+1]-B[c];
	}

    inline float& fxl() {return value_scaledf[0];}
    inline float& fyl() {return value_scaledf[1];}
    inline float& cxl() {return value_scaledf[2];}
    inline float& cyl() {return value_scaledf[3];}
};

struct FrameHessian
{
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    FrameShell* shell;
    float ab_exposure;

    Eigen::Vector3f* dI;				//!< 图像导数 // trace, fine tracking. Used for direction select (not for gradient histograms etc.)
	Eigen::Vector3f* dIp[PYR_LEVELS];	 // coarse tracking / coarse initializer. NAN in [0] only.
    float* absSquaredGrad[PYR_LEVELS];

	inline FrameHessian()
	{
		
	};

    void makeImages(float* color, CalibHessian* HCalib);
};


