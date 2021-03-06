#pragma once
#include "util/globalutil.h"
#include <iostream>
#include <fstream>
#include "FullSystem/HessianBlocks.h"
#include "IOWrapper/Output3DWrapper.h"
#include "boost/thread.hpp"
#include "util/ImageAndExposure.h"
#include "FullSystem/CoarseInitializer.h"
#include "OptimizationBackend/EnergyFunctional.h"
#include "ImmaturePoint.h"
#include "FullSystem/CoarseTracker.h"

using namespace dso;

class FullSystem
{
public:
	FullSystem();

	void setGammaFunction(float* BInv);
	void setPrecalcValues();
	void addActiveFrame( ImageAndExposure* image, int id );

	bool linearizeOperation;

	bool initFailed; 
	bool initialized;
	bool isLost;
	std::vector<IOWrap::Output3DWrapper*> outputWrapper;

protected:
	void initializeFromInitializer(FrameHessian* fh);
	void deliverTrackedFrame(FrameHessian* fh, bool needKF);
	void traceNewCoarse(FrameHessian* fh);

	EnergyFunctional* ef;			//!< 能量方程

	boost::mutex trackMutex;

	boost::mutex mapMutex;

	// mutex for camToWorl's in shells (these are always in a good configuration).
	boost::mutex shellPoseMutex;

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

	std::vector<FrameShell*> allKeyFramesHistory;
	std::vector<FrameHessian*> frameHessians;	//!< 关键帧 // ONLY changed in marginalizeFrame and addFrame.

	CoarseTracker* coarseTracker_forNewKF;			// set as as reference. protected by [coarseTrackerSwapMutex].
	CoarseTracker* coarseTracker;					// always used to track new frames. protected by [trackMutex].
	
	int lastRefStopID;

	void makeKeyFrame( FrameHessian* fh);
};