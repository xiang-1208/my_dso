#pragma once
#include <Eigen/Core>

class ImageAndExposure
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
  	float* image;			// irradiance. between 0 and 256
	int w,h;				// width and height;
	double timestamp;
	float exposure_time;	// exposure time in ms.

    inline ImageAndExposure(int w_, int h_, double timestamp_=0) : w(w_), h(h_), timestamp(timestamp_)
	{
		image = new float[w*h];
		exposure_time=1;
	}  

	inline void copyMetaTo(ImageAndExposure &other)
	{
		other.exposure_time = exposure_time;
	}

	inline ~ImageAndExposure()
	{
		delete[] image;
	} 
};
