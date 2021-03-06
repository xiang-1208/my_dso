#include "util/globalutil.h"

namespace dso
{
	int pyrLevelsUsed = PYR_LEVELS;

    float wM3G;
    float hM3G;
    int wG[PYR_LEVELS], hG[PYR_LEVELS];
    float fxG[PYR_LEVELS], fyG[PYR_LEVELS],
    	  cxG[PYR_LEVELS], cyG[PYR_LEVELS];
    float fxiG[PYR_LEVELS], fyiG[PYR_LEVELS],
    	  cxiG[PYR_LEVELS], cyiG[PYR_LEVELS];
    Eigen::Matrix3f KG[PYR_LEVELS], KiG[PYR_LEVELS];

    std::string vignette = "/home/xiang/dataset/tum/sequence_01/vignette.png";
    std::string gammaCalib = "/home/xiang/dataset/tum/sequence_01/pcalib.txt";
    std::string source = "/home/xiang/dataset/tum/sequence_01/images.zip";
    std::string calib = "/home/xiang/dataset/tum/sequence_01/camera.txt";
    // 2022.1.12
    std::string undistortion = "";
    // 2022.1.12
    std::string extrinsic = "";

	int mode=0;

    bool useSampleOutput=false;
    int start=0;
    int end=100000; 

    bool reverse = false;
    bool setting_logStuff = true;   

    // 0 = nothing.
    // 1 = apply inv. response.
    // 2 = apply inv. response & remove V.
    int setting_photometricCalibration = 2;
    float setting_affineOptModeA = 1e12; //-1: fix. >=0: optimize (with prior, if > 0).
    float setting_affineOptModeB = 1e8; //-1: fix. >=0: optimize (with prior, if > 0).
    float setting_minGradHistAdd = 7;
	float setting_minGradHistCut = 0.5;
	float setting_gradDownweightPerLevel = 0.75;
	bool  setting_selectDirectionDistribution = true;
	int setting_gammaWeightsPixelSelect = 1; // 1 = use original intensity for pixel selection; 0 = use gamma-corrected intensity.

	float setting_maxPixSearch = 0.027; // max length of the ep. line segment searched during immature point tracking. relative to image resolution.
	float setting_trace_slackInterval = 1.5;			// if pixel-interval is smaller than this, leave it be.
	float setting_trace_stepsize = 1.0;				// stepsize for initial discrete search.
	float setting_trace_minImprovementFactor = 2;		// if pixel-interval is smaller than this, leave it be.
	int setting_minTraceTestRadius = 2;
	int setting_trace_GNIterations = 3;				// max # GN iterations
	float setting_trace_GNThreshold = 0.1;				// GN stop after this stepsize.
	float setting_trace_extraSlackOnTH = 1.2;			// for energy-based outlier check, be slightly more relaxed by this factor.

	float setting_idepthFixPrior = 50*50;
	float setting_idepthFixPriorMargFac = 600*600;
	float setting_initialRotPrior = 1e11;
	float setting_initialTransPrior = 1e10;
	float setting_initialAffBPrior = 1e14;
	float setting_initialAffAPrior = 1e14;
	float setting_initialCalibHessian = 5e9;

	/* some modes for solving the resulting linear system (e.g. orthogonalize wrt. unobservable dimensions) */
	int setting_solverMode = SOLVER_FIX_LAMBDA | SOLVER_ORTHOGONALIZE_X_LATER;  // ??????????????????lambda

    float playbackSpeed=0;	// 0 for linearize (play as fast as possible, while sequentializing tracking & mapping). otherwise, factor on timestamps.
    bool disableAllDisplay = false;

    bool setting_render_display3D = true;
	bool setting_debugout_runquiet = false;

    float setting_desiredImmatureDensity = 1500; // immature points per frame
    float setting_desiredPointDensity = 2000; // aimed total points in the active window.
    float setting_minPointsRemaining = 0.05;  // marg a frame if less than X% points remain.
    float setting_maxLogAffFacInWindow = 0.7; // marg a frame if factor between intensities to current frame is larger than 1/X or X.

    int   setting_minFrames = 5; // min frames in window.
    int   setting_maxFrames = 7; // max frames in window.

    float setting_kfGlobalWeight = 1;   // general weight on threshold, the larger the more KF's are taken (e.g., 2 = double the amount of KF's).
    float setting_maxAffineWeight= 2;

	float setting_outlierTH = 12*12;					// higher -> less strict
	float setting_outlierTHSumComponent = 50*50; 		// higher -> less strong gradient-based reweighting .
	float setting_overallEnergyTHWeight = 1;

	bool setting_render_displayResidual = true;
    bool setting_render_displayVideo = true;
    bool setting_render_displayDepth = true;
	bool setting_render_displayCoarseTrackingFull=false;
	bool setting_render_renderWindowFrames=true;
	bool setting_render_plotTrackingFull = false;
	bool setting_fullResetRequested = false;

	bool preload=false;

	float setting_huberTH = 9; // Huber Threshold

	int sparsityFactor = 5;	// not actually a setting, only some legacy stuff for coarse initializer.

	float benchmark_varNoise = 0;
	float benchmark_varBlurNoise = 0;
	int benchmark_noiseGridsize = 3;

	float freeDebugParam1 = 1;
	float freeDebugParam2 = 1;
	float freeDebugParam3 = 1;
	float freeDebugParam4 = 1;
	float freeDebugParam5 = 1;

	bool goStepByStep = false;

	void handleKey(char k)
	{
		char kkk = k;
		switch(kkk)
		{
		case 'd': case 'D':
			freeDebugParam5 = ((int)(freeDebugParam5+1))%10;
			printf("new freeDebugParam5: %f!\n", freeDebugParam5);
			break;
		case 's': case 'S':
			freeDebugParam5 = ((int)(freeDebugParam5-1+10))%10;
			printf("new freeDebugParam5: %f!\n", freeDebugParam5);
			break;
		}
	
	}


int staticPattern[10][40][2] = {
		{{0,0}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// .
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{0,-1},	  {-1,0},	   {0,0},	    {1,0},	     {0,1}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// +
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{-1,-1},	  {1,1},	   {0,0},	    {-1,1},	     {1,-1}, 	  {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},	// x
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{-1,-1},	  {-1,0},	   {-1,1},		{-1,0},		 {0,0},		  {0,1},	   {1,-1},		{1,0},		 {1,1},       {-100,-100},	// full-tight
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-100,-100},	// full-spread-9
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-2,-2},   // full-spread-13
		 {-2,2},      {2,-2},      {2,2},       {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{-2,-2},     {-2,-1}, {-2,-0}, {-2,1}, {-2,2}, {-1,-2}, {-1,-1}, {-1,-0}, {-1,1}, {-1,2}, 										// full-25
		 {-0,-2},     {-0,-1}, {-0,-0}, {-0,1}, {-0,2}, {+1,-2}, {+1,-1}, {+1,-0}, {+1,1}, {+1,2},
		 {+2,-2}, 	  {+2,-1}, {+2,-0}, {+2,1}, {+2,2}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{1,1},		 {0,2},       {-2,-2},   // full-spread-21
		 {-2,2},      {2,-2},      {2,2},       {-3,-1},     {-3,1},      {3,-1}, 	   {3,1},       {1,-3},      {-1,-3},     {1,3},
		 {-1,3},      {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{0,-2},	  {-1,-1},	   {1,-1},		{-2,0},		 {0,0},		  {2,0},	   {-1,1},		{0,2},		 {-100,-100}, {-100,-100},	// 8 for SSE efficiency
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100},
		 {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}, {-100,-100}},

		{{-4,-4},     {-4,-2}, {-4,-0}, {-4,2}, {-4,4}, {-2,-4}, {-2,-2}, {-2,-0}, {-2,2}, {-2,4}, 										// full-45-SPREAD
		 {-0,-4},     {-0,-2}, {-0,-0}, {-0,2}, {-0,4}, {+2,-4}, {+2,-2}, {+2,-0}, {+2,2}, {+2,4},
		 {+4,-4}, 	  {+4,-2}, {+4,-0}, {+4,2}, {+4,4}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200},
		 {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}, {-200,-200}},
};





};