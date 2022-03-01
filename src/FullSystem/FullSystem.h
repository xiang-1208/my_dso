#pragma once
#include "util/globalutil.h"
#include <iostream>
#include <fstream>
#include "FullSystem/HessianBlocks.h"
#include "IOWrapper/Output3DWrapper.h"
#include "boost/thread.hpp"
#include "util/ImageAndExposure.h"
#include "FullSystem/CoarseInitializer.h"

using namespace dso;

class FullSystem
{
public:
	FullSystem();

	void setGammaFunction(float* BInv);
	void addActiveFrame( ImageAndExposure* image, int id );

	bool linearizeOperation;
	bool initialized;
	bool isLost;
	std::vector<IOWrap::Output3DWrapper*> outputWrapper;

protected:
	boost::mutex trackMutex;

	CalibHessian Hcalib;

	std::ofstream* calibLog;
	std::ofstream* numsLog;
	std::ofstream* errorsLog;
	std::ofstream* eigenAllLog;
	std::ofstream* eigenPLog;
	std::ofstream* eigenALog;
	std::ofstream* DiagonalLog;
	std::ofstream* variancesLog;
	std::ofstream* nullspacesLog;

	std::ofstream* coarseTrackingLog;

	float* selectionMap;

	std::vector<FrameShell*> allFrameHistory;
	CoarseInitializer* coarseInitializer;
	
};