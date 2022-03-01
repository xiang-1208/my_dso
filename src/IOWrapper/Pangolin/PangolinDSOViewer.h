#pragma once
#include <Eigen/Core>
#include <util/globalutil.h>
#include <pangolin/pangolin.h>
#include "boost/thread.hpp"
#include "IOWrapper/Output3DWrapper.h"
#include "util/MinimalImage.h"
#include "KeyFrameDisplay.h"
#include <deque>

using namespace dso;

namespace IOWrap
{

struct GraphConnection
{
	KeyFrameDisplay* from;
	KeyFrameDisplay* to;
	int fwdMarg, bwdMarg, fwdAct, bwdAct;
};

class PangolinDSOViewer : public Output3DWrapper
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    PangolinDSOViewer(int w, int h, bool startRunThread=true);
    void run();

    virtual void pushLiveFrame(FrameHessian* image) override;

    virtual void join() override;
private:

    void drawConstraints();
    void reset_internal();

    int w;
    int h;
    boost::mutex openImagesMutex;
    boost::mutex model3DMutex;

	float settings_scaledVarTH;
	float settings_absVarTH;
	int settings_pointCloudMode;
	float settings_minRelBS;
	int settings_sparsity;

    boost::thread runThread;

    MinimalImageB3* internalVideoImg;
	MinimalImageB3* internalKFImg;
	MinimalImageB3* internalResImg;

    KeyFrameDisplay* currentCam;
    std::vector<KeyFrameDisplay*> keyframes;
    std::vector<Vec3f,Eigen::aligned_allocator<Vec3f>> allFramePoses;
    std::map<int, KeyFrameDisplay*> keyframesByKFID;
    std::vector<GraphConnection,Eigen::aligned_allocator<GraphConnection>> connections;
    std::deque<float> lastNTrackingMs;
	std::deque<float> lastNMappingMs;

    bool videoImgChanged,kfImgChanged,resImgChanged;
    bool needReset;
    bool running;
    bool settings_showKFCameras;
    bool settings_showCurrentCamera;
	bool settings_showTrajectory;
	bool settings_showFullTrajectory;
	bool settings_showActiveConstraints;
	bool settings_showAllConstraints;

};

}