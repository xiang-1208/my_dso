#pragma once
#include "Eigen/Core"
#include <util/NumType.h>

template<typename T>
class MinimalImage
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	int w;
	int h;
	T* data;

	/*
	 * creates minimal image with own memory
	 */
	inline MinimalImage(int w_, int h_) : w(w_), h(h_)
	{
		data = new T[w*h];
		ownData=true;
	}

	/*
	 * creates minimal image wrapping around existing memory
	 */
	inline MinimalImage(int w_, int h_, T* data_) : w(w_), h(h_)
	{
		data = data_;
		ownData=false;
	}

	inline void setBlack()
	{
		memset(data, 0, sizeof(T)*w*h);
	}

	inline T at(int i) {return data[i];}
	inline T& at(int x, int y) {return data[(int)x+((int)y)*w];}

	inline void setPixelCirc(const int &u, const int &v, T val)
	{
		for(int i=-3;i<=3;i++)
		{
			at(u+3,v+i) = val;
			at(u-3,v+i) = val;
			at(u+2,v+i) = val;
			at(u-2,v+i) = val;

			at(u+i,v-3) = val;
			at(u+i,v+3) = val;
			at(u+i,v-2) = val;
			at(u+i,v+2) = val;
		}
	}
    
private:
	bool ownData;
};

typedef MinimalImage<unsigned char> MinimalImageB;
typedef MinimalImage<Vec3b> MinimalImageB3;
typedef MinimalImage<Vec3f> MinimalImageF3;