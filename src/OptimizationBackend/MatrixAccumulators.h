#pragma once
#include "util/NumType.h"

#if !defined(__SSE3__) && !defined(__SSE2__) && !defined(__SSE1__)
#include "SSE2NEON.h"
#endif


template<int i>
class AccumulatorX
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

  	Eigen::Matrix<float,i,1> A;
  	Eigen::Matrix<float,i,1> A1k;
  	Eigen::Matrix<float,i,1> A1m;	
	size_t num;

  	inline void initialize()
  	{
  	  	A.setZero();
  	  	A1k.setZero();
  	  	A1m.setZero();
  	  	num = numIn1 = numIn1k = numIn1m = 0;
  	}

	inline void finish()
  	{
		shiftUp(true);
		num = numIn1+numIn1k+numIn1m;
  	}


  	inline void update(const Eigen::Matrix<float,i,1> &L, float w)
  	{
		A += w*L;
		numIn1++;
		shiftUp(false);
  	}

  	inline void updateNoWeight(const Eigen::Matrix<float,i,1> &L)
  	{
		A += L;
		numIn1++;
		shiftUp(false);
  	}
	
private:
	float numIn1, numIn1k, numIn1m;

  	void shiftUp(bool force)
  	{
		if(numIn1 > 1000 || force)
		{
		  	A1k += A;
		  	A.setZero();
		  	numIn1k+=numIn1;
		  	numIn1=0;
		}
		if(numIn1k > 1000 || force)
		{
		  	A1m += A1k;
		  	A1k.setZero();
		  	numIn1m+=numIn1k;
		  	numIn1k=0;
		}
  	}
};

class Accumulator11
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

    float A;
    size_t num;

    inline void initialize()
    {
      A=0;
      memset(SSEData,0, sizeof(float)*4*1);
      memset(SSEData1k,0, sizeof(float)*4*1);
      memset(SSEData1m,0, sizeof(float)*4*1);
      num = numIn1 = numIn1k = numIn1m = 0;
    }   

    inline void updateSingle(const float val)
    {
        SSEData[0] += val;
        num++; numIn1++;
        shiftUp(false);
    }

    inline void finish()
    {
	    shiftUp(true); // 都进位到 m
	    A=SSEData1m[0+0] + SSEData1m[0+1] + SSEData1m[0+2] + SSEData1m[0+3];
    }

private:
    EIGEN_ALIGN16 float SSEData[4*1];  // 16字节
    EIGEN_ALIGN16 float SSEData1k[4*1];
    EIGEN_ALIGN16 float SSEData1m[4*1];
    float numIn1, numIn1k, numIn1m; 

    void shiftUp(bool force)
    {
        // 大于1000, 相加则进位到 k 
        if(numIn1 > 1000 || force)
        {
            _mm_store_ps(SSEData1k, _mm_add_ps(_mm_load_ps(SSEData),_mm_load_ps(SSEData1k)));
            numIn1k+=numIn1; numIn1=0;
            memset(SSEData,0, sizeof(float)*4*1);
        }

	    if(numIn1k > 1000 || force)
	    {
	    	_mm_store_ps(SSEData1m, _mm_add_ps(_mm_load_ps(SSEData1k),_mm_load_ps(SSEData1m)));
	    	numIn1m+=numIn1k;  numIn1k=0;
	    	memset(SSEData1k,0, sizeof(float)*4*1);
	    }
    }
};

class Accumulator9
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;   

    Mat99f H;
    Vec9f b;
    size_t num;

    inline void initialize()
    {
    	H.setZero();
    	b.setZero(); 
        memset(SSEData,0, sizeof(float)*4*45);  // 会对128位, 16字节进行对齐, 因此每个数用4个float存
        memset(SSEData1k,0, sizeof(float)*4*45);
        memset(SSEData1m,0, sizeof(float)*4*45);
        num = numIn1 = numIn1k = numIn1m = 0;     
    }   

    // 计算一个9维向量相乘, 得到9*9矩阵
    inline void updateSSE(
		const __m128 J0,const __m128 J1,
		const __m128 J2,const __m128 J3,
		const __m128 J4,const __m128 J5,
		const __m128 J6,const __m128 J7,
		const __m128 J8)      
    {
 		// 一共45个值
	    float* pt=SSEData;
		// 第一行9个值
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J0))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J1))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J2))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J3))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J4))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J0,J8))); pt+=4;
	    	// 第二行8个, 因为对称
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J1))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J2))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J3))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J4))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J1,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J2))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J3))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J4))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J2,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J3))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J4))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J3,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J4,J4))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J4,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J4,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J4,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J4,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J5,J5))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J5,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J5,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J5,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J6,J6))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J6,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J6,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J7,J7))); pt+=4;
	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J7,J8))); pt+=4;

	     _mm_store_ps(pt, _mm_add_ps(_mm_load_ps(pt),_mm_mul_ps(J8,J8))); pt+=4;

	     num+=4;
	     numIn1++;  // 乘一次加一
	     shiftUp(false);       
    }  

    // 不使用_m128来计算
    inline void updateSingle(
		const float J0,const float J1,
		const float J2,const float J3,
		const float J4,const float J5,
		const float J6,const float J7,
		const float J8, int off=0)
    {
	    float* pt=SSEData+off;
	    *pt += J0*J0; pt+=4;
	    *pt += J1*J0; pt+=4;
	    *pt += J2*J0; pt+=4;
	    *pt += J3*J0; pt+=4;
	    *pt += J4*J0; pt+=4;
	    *pt += J5*J0; pt+=4;
	    *pt += J6*J0; pt+=4;
	    *pt += J7*J0; pt+=4;
	    *pt += J8*J0; pt+=4;


	    *pt += J1*J1; pt+=4;
	    *pt += J2*J1; pt+=4;
	    *pt += J3*J1; pt+=4;
	    *pt += J4*J1; pt+=4;
	    *pt += J5*J1; pt+=4;
	    *pt += J6*J1; pt+=4;
	    *pt += J7*J1; pt+=4;
	    *pt += J8*J1; pt+=4;


	    *pt += J2*J2; pt+=4;
	    *pt += J3*J2; pt+=4;
	    *pt += J4*J2; pt+=4;
	    *pt += J5*J2; pt+=4;
	    *pt += J6*J2; pt+=4;
	    *pt += J7*J2; pt+=4;
	    *pt += J8*J2; pt+=4;


	    *pt += J3*J3; pt+=4;
	    *pt += J4*J3; pt+=4;
	    *pt += J5*J3; pt+=4;
	    *pt += J6*J3; pt+=4;
	    *pt += J7*J3; pt+=4;
	    *pt += J8*J3; pt+=4;


	    *pt += J4*J4; pt+=4;
	    *pt += J5*J4; pt+=4;
	    *pt += J6*J4; pt+=4;
	    *pt += J7*J4; pt+=4;
	    *pt += J8*J4; pt+=4;

	    *pt += J5*J5; pt+=4;
	    *pt += J6*J5; pt+=4;
	    *pt += J7*J5; pt+=4;
	    *pt += J8*J5; pt+=4;


	    *pt += J6*J6; pt+=4;
	    *pt += J7*J6; pt+=4;
	    *pt += J8*J6; pt+=4;


	    *pt += J7*J7; pt+=4;
	    *pt += J8*J7; pt+=4;

	    *pt += J8*J8; pt+=4;

	    num++;
	    numIn1++;
	    shiftUp(false);
    }

	// 不使用对齐加速的, 带有权重的
  	inline void updateSingleWeighted(
		float J0, float J1,
		float J2, float J3,
		float J4, float J5,
		float J6, float J7,
		float J8, float w,
		int off=0)
  	{
	  	float* pt=SSEData+off;
	  	*pt += J0*J0*w; pt+=4; J0*=w;
	  	*pt += J1*J0; pt+=4;
	  	*pt += J2*J0; pt+=4;
	  	*pt += J3*J0; pt+=4;
	  	*pt += J4*J0; pt+=4;
	  	*pt += J5*J0; pt+=4;
	  	*pt += J6*J0; pt+=4;
	  	*pt += J7*J0; pt+=4;
	  	*pt += J8*J0; pt+=4;
	
	
	  	*pt += J1*J1*w; pt+=4; J1*=w;
	  	*pt += J2*J1; pt+=4;
	  	*pt += J3*J1; pt+=4;
	  	*pt += J4*J1; pt+=4;
	  	*pt += J5*J1; pt+=4;
	  	*pt += J6*J1; pt+=4;
	  	*pt += J7*J1; pt+=4;
	  	*pt += J8*J1; pt+=4;
	
	
	  	*pt += J2*J2*w; pt+=4; J2*=w;
	  	*pt += J3*J2; pt+=4;
	  	*pt += J4*J2; pt+=4;
	  	*pt += J5*J2; pt+=4;
	  	*pt += J6*J2; pt+=4;
	  	*pt += J7*J2; pt+=4;
	  	*pt += J8*J2; pt+=4;
	
	
	  	*pt += J3*J3*w; pt+=4; J3*=w;
	  	*pt += J4*J3; pt+=4;
	  	*pt += J5*J3; pt+=4;
	  	*pt += J6*J3; pt+=4;
	  	*pt += J7*J3; pt+=4;
	  	*pt += J8*J3; pt+=4;
	
	
	  	*pt += J4*J4*w; pt+=4; J4*=w;
	  	*pt += J5*J4; pt+=4;
	  	*pt += J6*J4; pt+=4;
	  	*pt += J7*J4; pt+=4;
	  	*pt += J8*J4; pt+=4;
	
	  	*pt += J5*J5*w; pt+=4; J5*=w;
	  	*pt += J6*J5; pt+=4;
	  	*pt += J7*J5; pt+=4;
	  	*pt += J8*J5; pt+=4;
	
	
	  	*pt += J6*J6*w; pt+=4; J6*=w;
	  	*pt += J7*J6; pt+=4;
	  	*pt += J8*J6; pt+=4;
	
	
	  	*pt += J7*J7*w; pt+=4; J7*=w;
	  	*pt += J8*J7; pt+=4;
	
	  	*pt += J8*J8*w; pt+=4;
	
	  	num++;
	  	numIn1++;
	  	shiftUp(false);
  	}


    inline void finish()
    {
	    H.setZero();
	    shiftUp(true);  // 强制进位到m
	    assert(numIn1==0);
	    assert(numIn1k==0);

	    int idx=0;
	    //* H矩阵是对称的, 只有45个数值
	    for(int r=0;r<9;r++)
	    	for(int c=r;c<9;c++)
	    	{
	    		float d = SSEData1m[idx+0] + SSEData1m[idx+1] + SSEData1m[idx+2] + SSEData1m[idx+3];
	    		H(r,c) = H(c,r) = d;
	    		idx+=4;
	    	}
	      assert(idx==4*45);
    }

private:
    EIGEN_ALIGN16 float SSEData[4*45];
    EIGEN_ALIGN16 float SSEData1k[4*45];
    EIGEN_ALIGN16 float SSEData1m[4*45];
    float numIn1, numIn1k, numIn1m;

    void shiftUp(bool force)
    {
	    if(numIn1 > 1000 || force)
	    {
		    for(int i=0;i<45;i++)
		    	_mm_store_ps(SSEData1k+4*i, _mm_add_ps(_mm_load_ps(SSEData+4*i),_mm_load_ps(SSEData1k+4*i)));
		    numIn1k+=numIn1;
		    numIn1=0;
		    memset(SSEData,0, sizeof(float)*4*45);
	    }

	    if(numIn1k > 1000 || force)
	    {
		    for(int i=0;i<45;i++)
		    	_mm_store_ps(SSEData1m+4*i, _mm_add_ps(_mm_load_ps(SSEData1k+4*i),_mm_load_ps(SSEData1m+4*i)));
		    numIn1m+=numIn1k;
		    numIn1k=0;
		    memset(SSEData1k,0, sizeof(float)*4*45);
	    }
    }
};