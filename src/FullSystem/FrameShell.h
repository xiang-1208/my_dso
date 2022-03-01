#pragma once
#include <Eigen/Core>
#include "util/NumType.h"

class FrameShell
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    int id; 	// INTERNAL ID, starting at zero.
    int incoming_id;
    double timestamp;

    SE3 camToTrackingRef;
    FrameShell* trackingRef;

    bool poseValid;
    AffLight aff_g2l;       //
    SE3 camToWorld;

    // statisitcs
    int marginalizedAt;     //边缘化
    double movedByOpt;
	int statistics_outlierResOnThis;
	int statistics_goodResOnThis;


    inline FrameShell()
    {
        id=0;
        poseValid=true;
        camToWorld = SE3();
        timestamp=0;
        marginalizedAt=-1;
        movedByOpt=0;
        statistics_outlierResOnThis=statistics_goodResOnThis=0;
		trackingRef=0;
		camToTrackingRef = SE3();        
    }
};