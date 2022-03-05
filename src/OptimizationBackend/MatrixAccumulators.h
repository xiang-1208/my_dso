#pragma once
#include "util/NumType.h"

#if !defined(__SSE3__) && !defined(__SSE2__) && !defined(__SSE1__)
#include "SSE2NEON.h"
#endif

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

private:
    EIGEN_ALIGN16 float SSEData[4*45];
    EIGEN_ALIGN16 float SSEData1k[4*45];
    EIGEN_ALIGN16 float SSEData1m[4*45];
    float numIn1, numIn1k, numIn1m;
};