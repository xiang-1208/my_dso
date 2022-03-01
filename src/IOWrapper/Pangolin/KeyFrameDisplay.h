#pragma once
#include <util/NumType.h>
#include <util/globalutil.h>
#include <pangolin/pangolin.h>
#include <Eigen/Core>

using namespace dso;

template<int ppp>
struct InputPointSparse
{
	float u;
	float v;
	float idpeth;
	float idepth_hessian;
	float relObsBaseline;
	int numGoodRes;
	unsigned char color[ppp];
	unsigned char status;
};

class KeyFrameDisplay
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	KeyFrameDisplay();

	void drawCam(float lineWidth, float* color, float sizeFactor);
	bool refreshPC(bool canRefresh, float scaledTH, float absTH, int mode, float minBS, int sparsity);
	void drawPC(float pointSize);

	int id;
	bool active;
	SE3 camToWorld;

protected:
	int width, height;
	float fx,fy,cx,cy;
	float fxi,fyi,cxi,cyi;

	int numSparsePoints;
	int numSparseBufferSize;
    InputPointSparse<MAX_RES_PER_POINT>* originalInputSparse;

	float my_scaledTH, my_absTH, my_scale;
	int my_sparsifyFactor;
	int my_displayMode;
	float my_minRelBS;
	bool needRefresh;
	
	bool bufferValid;
	int numGLBufferPoints;
	int numGLBufferGoodPoints;

	pangolin::GlBuffer vertexBuffer;
	pangolin::GlBuffer colorBuffer;
};