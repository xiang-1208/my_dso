#pragma once
#include <string>
#include <Eigen/Dense>
#include <iostream>

#define PYR_LEVELS 6

namespace dso
{
    #define SOLVER_SVD (int)1
    #define SOLVER_ORTHOGONALIZE_SYSTEM (int)2
    #define SOLVER_ORTHOGONALIZE_POINTMARG (int)4
    #define SOLVER_ORTHOGONALIZE_FULL (int)8
    #define SOLVER_SVD_CUT7 (int)16
    #define SOLVER_REMOVE_POSEPRIOR (int)32
    #define SOLVER_USE_GN (int)64
    #define SOLVER_FIX_LAMBDA (int)128
    #define SOLVER_ORTHOGONALIZE_X (int)256
    #define SOLVER_MOMENTUM (int)512
    #define SOLVER_STEPMOMENTUM (int)1024
    #define SOLVER_ORTHOGONALIZE_X_LATER (int)2048

    #define patternNum 8

    extern int pyrLevelsUsed;
    
    
    extern float wM3G;
    extern float hM3G;
    extern int wG[PYR_LEVELS], hG[PYR_LEVELS];
    extern float fxG[PYR_LEVELS], fyG[PYR_LEVELS],
	        cxG[PYR_LEVELS], cyG[PYR_LEVELS];
    extern float fxiG[PYR_LEVELS], fyiG[PYR_LEVELS],
            cxiG[PYR_LEVELS], cyiG[PYR_LEVELS];
    extern Eigen::Matrix3f KG[PYR_LEVELS], KiG[PYR_LEVELS];

    extern std::string vignette;
    extern std::string gammaCalib;
    extern std::string source;
    extern std::string calib;
 
    extern std::string undistortion;

    extern std::string extrinsic;
    extern bool useSampleOutput;
    extern int mode;
    extern int start;
    extern int end; 
 
    extern bool reverse;
    extern bool setting_logStuff; 

    extern float setting_maxPixSearch;
    extern float setting_trace_slackInterval;
    extern float setting_trace_stepsize;
    extern float setting_trace_minImprovementFactor;
    extern int setting_minTraceTestRadius;
    extern int setting_trace_GNIterations;
    extern float setting_trace_GNThreshold;

	extern float setting_idepthFixPrior;
	extern float setting_idepthFixPriorMargFac;
	extern float setting_initialRotPrior;
	extern float setting_initialTransPrior;
	extern float setting_initialAffBPrior;
	extern float setting_initialAffAPrior;
	extern float setting_initialCalibHessian; 
    extern float setting_trace_extraSlackOnTH; 


    extern int setting_photometricCalibration ;
    extern float setting_affineOptModeA; //-1: fix. >=0: optimize (with prior, if > 0).
    extern float setting_affineOptModeB; //-1: fix. >=0: optimize (with prior, if > 0).
    extern float setting_minGradHistAdd;
    extern int setting_gammaWeightsPixelSelect;

    extern float playbackSpeed;
    extern bool disableAllDisplay;

    extern float setting_outlierTH;

    extern float setting_huberTH;

    extern float setting_minGradHistAdd;
    extern float setting_minGradHistCut;
    extern float setting_gradDownweightPerLevel;
    extern bool  setting_selectDirectionDistribution;

    extern bool setting_render_display3D;
    extern float setting_desiredImmatureDensity;
    extern float setting_desiredPointDensity;
    extern float setting_minPointsRemaining;
    extern float setting_maxLogAffFacInWindow;

    extern int setting_solverMode;

    extern int setting_minFrames;
    extern int setting_maxFrames;

    extern float setting_kfGlobalWeight;
    extern float setting_maxAffineWeight;
    
    extern int staticPattern[10][40][2];

    #define patternNum 8
    #define patternP staticPattern[8]   
    #define patternPadding 2

    extern bool setting_render_displayResidual;
    extern bool setting_render_displayVideo;
    extern bool setting_render_displayDepth;
    extern bool setting_render_displayCoarseTrackingFull;
    extern bool setting_render_renderWindowFrames;
    extern bool setting_render_plotTrackingFull ;
    extern bool setting_fullResetRequested;

    extern bool setting_debugout_runquiet;
    extern float setting_outlierTHSumComponent;

    extern float setting_overallEnergyTHWeight;

    extern bool preload;

    extern int sparsityFactor;

    extern float benchmark_varNoise;
    extern float benchmark_varBlurNoise;
    extern int benchmark_noiseGridsize;

	extern float freeDebugParam1;
	extern float freeDebugParam2;
	extern float freeDebugParam3;
	extern float freeDebugParam4;
	extern float freeDebugParam5;

    extern bool goStepByStep;

    void handleKey(char k);
    
}