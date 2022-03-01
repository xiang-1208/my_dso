#pragma once
#include "util/NumType.h"
#include "util/globalutil.h"
#include "FullSystem/HessianBlocks.h"
#include "util/MinimalImage.h"
#include "IOWrapper/ImageDisplay.h"

using namespace dso;

enum PixelSelectorStatus {PIXSEL_VOID=0, PIXSEL_1, PIXSEL_2, PIXSEL_3};

class PixelSelector
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    PixelSelector(int w, int h);
	int makeMaps(
			const FrameHessian* const fh,
			float* map_out, float density, int recursionsLeft=1, bool plot=false, float thFactor=1);
    

    int currentPotential;

private:
    unsigned char* randomPattern;
    const FrameHessian* gradHistFrame;
    int thsStep;
    int* gradHist;
    float* ths;
    float* thsSmoothed;

    void makeHists(const FrameHessian* const fh);
    Eigen::Vector3i select(const FrameHessian* const fh,
		float* map_out, int pot, float thFactor);
};