#include "FullSystem.h"

int CalibHessian::instanceCounter=0;

FullSystem::FullSystem()
{
    int retstat=0;
    if (setting_logStuff)
    {
		retstat += system("rm -rf logs");
		retstat += system("mkdir logs");

		retstat += system("rm -rf mats");
		retstat += system("mkdir mats");   

        calibLog = new std::ofstream;
        calibLog->open("logs/caliblog.txt",std::ios::trunc | std::ios::out);
        calibLog->precision(12);

        numsLog = new std::ofstream;
        numsLog->open("logs/numsLog.txt",std::ios::trunc | std::ios::out);
        numsLog->precision(10);

        errorsLog = new std::ofstream;
        errorsLog->open("logs/errorsLog.txt",std::ios::trunc | std::ios::out);
        errorsLog->precision(10);

        eigenAllLog = new std::ofstream;
        eigenAllLog->open("logs/eigenAllLog.txt",std::ios::trunc | std::ios::out);
        eigenAllLog->precision(10);

        eigenPLog = new std::ofstream;
        eigenPLog->open("logs/eigenPLog.txt",std::ios::trunc | std::ios::out);
        eigenPLog->precision(10);

        eigenALog = new std::ofstream;
        eigenALog->open("logs/eigenALog.txt",std::ios::trunc | std::ios::out);
        eigenALog->precision(10);   

        DiagonalLog = new std::ofstream;
        DiagonalLog->open("logs/DiagonalLog.txt",std::ios::trunc | std::ios::out);
        DiagonalLog->precision(10);

        variancesLog = new std::ofstream;
        variancesLog->open("logs/variancesLog.txt",std::ios::trunc | std::ios::out);
        variancesLog->precision(10);

        nullspacesLog = new std::ofstream;
        nullspacesLog->open("logs/nullspacesLog.txt",std::ios::trunc | std::ios::out);
        nullspacesLog->precision(10);     
    }
    else
    {
	    calibLog = 0;
	    numsLog = 0;
	    errorsLog = 0;
	    eigenAllLog = 0;
	    eigenPLog = 0;
	    eigenALog = 0;
	    DiagonalLog = 0;
	    variancesLog = 0;
	    nullspacesLog = 0;       
    }

    assert(retstat!=293847);

    coarseInitializer = new CoarseInitializer(wG[0], hG[0]);

    ef = new EnergyFunctional();

    //selectionMap = new float[wG[0]*hG[0]];
    initialized=false;
    isLost=false;
}

void FullSystem::setGammaFunction(float* BInv)
{
    if (BInv == 0)
        return;

    memcpy(Hcalib.Binv,BInv,sizeof(float)*256);

    for(int i=1;i<255;i++)
    {
        for (int s=1;s<255;s++)
        {
            if(BInv[s]<i && BInv[s+1]>i)
            {
				Hcalib.B[i] = s+(i - BInv[s]) / (BInv[s+1]-BInv[s]);
				break;
            }
        }
    }
    Hcalib.B[0]=0;
    Hcalib.B[255]=255;
}

void FullSystem::addActiveFrame( ImageAndExposure* image, int id )
{
    if(isLost) return;
    boost::unique_lock<boost::mutex> lock(trackMutex);

    FrameHessian* fh = new FrameHessian();
    FrameShell* shell = new FrameShell();
    shell->camToWorld = SE3();
    shell->aff_g2l = AffLight(0,0);

    shell->id = allFrameHistory.size();
    shell->marginalizedAt = shell->id;
    shell->timestamp = image->timestamp;
    shell->incoming_id = id;
	fh->shell = shell;     //!< 帧的"壳", 保存一些不变的,要留下来的量
	allFrameHistory.push_back(shell);  // 只把简略的shell存起来

	// 得到曝光时间
	fh->ab_exposure = image->exposure_time;
	// 构建金字塔以及计算梯度
    fh->makeImages(image->image, &Hcalib);

    //std::cout << "!!!!!!!!!!!!!!!!" <<std::endl;

    if(!initialized)
    {
        //初始化！！！初始化！！！
        if(coarseInitializer->frameID<0)    // first frame set. fh is kept by coarseInitializer.
        {
			//xiang finish 2022.2.14
            coarseInitializer->setFirst(&Hcalib, fh);
        }
        else if (coarseInitializer->trackFrame(fh, outputWrapper)) // if SNAPPED
        {
            std::cout << "初始化完成" <<std::endl;
            initializeFromInitializer(fh);
            lock.unlock();
        }
    }
    else
    {
    }
}

//@ 从初始化中提取出信息, 用于跟踪.
void FullSystem::initializeFromInitializer(FrameHessian* newFrame)
{
    boost::unique_lock<boost::mutex> lock(mapMutex);

    //[ ***step 1*** ] 把第一帧设置成关键帧, 加入队列, 加入EnergyFunctional
    FrameHessian* firstFrame = coarseInitializer->firstFrame;  // 第一帧增加进地图
	firstFrame->idx = frameHessians.size(); // 赋值给它id (0开始)
	frameHessians.push_back(firstFrame);  	// 地图内关键帧容器    
	firstFrame->frameID = allKeyFramesHistory.size();  	// 所有历史关键帧id
	allKeyFramesHistory.push_back(firstFrame->shell); 	// 所有历史关键帧
    ef->insertFrame(firstFrame, &Hcalib);
    setPrecalcValues();   		// 设置相对位姿预计算值

	//int numPointsTotal = makePixelStatus(firstFrame->dI, selectionMap, wG[0], hG[0], setting_desiredDensity);
	//int numPointsTotal = pixelSelector->makeMaps(firstFrame->dIp, selectionMap,setting_desiredDensity);

    firstFrame->pointHessians.reserve(wG[0]*hG[0]*0.2f); // 20%的点数目
	firstFrame->pointHessiansMarginalized.reserve(wG[0]*hG[0]*0.2f); // 被边缘化
	firstFrame->pointHessiansOut.reserve(wG[0]*hG[0]*0.2f); // 丢掉的点

    //[ ***step 2*** ] 求出平均尺度因子
    float sumID=1e-5, numID=1e-5;
	for(int i=0;i<coarseInitializer->numPoints[0];i++)
	{
		//? iR的值到底是啥
		sumID += coarseInitializer->points[0][i].iR; // 第0层点的中位值, 相当于
		numID++;
	}
    float rescaleFactor = 1 / (sumID / numID);  // 求出尺度因子

	// randomly sub-select the points I need.
	// 目标点数 / 实际提取点数    
    float keepPercentage = setting_desiredPointDensity / coarseInitializer->numPoints[0];

    if(!setting_debugout_runquiet)
        printf("Initialization: keep %.1f%% (need %d, have %d)!\n", 100*keepPercentage,
            (int)(setting_desiredPointDensity), coarseInitializer->numPoints[0] );

    //[ ***step 3*** ] 创建PointHessian, 点加入关键帧, 加入EnergyFunctional
    for(int i=0;i<coarseInitializer->numPoints[0];i++)
    {
        if(rand()/(float)RAND_MAX > keepPercentage) continue; // 如果提取的点比较少, 不执行; 提取的多, 则随机干掉

		Pnt* point = coarseInitializer->points[0]+i;
		ImmaturePoint* pt = new ImmaturePoint(point->u+0.5f,point->v+0.5f,firstFrame,point->my_type, &Hcalib);      

        if(!std::isfinite(pt->energyTH)) { delete pt; continue; }  // 点值无穷大  

        // 创建ImmaturePoint就为了创建PointHessian? 是为了接口统一吧
        pt->idepth_max=pt->idepth_min=1;
		PointHessian* ph = new PointHessian(pt, &Hcalib);
		delete pt;
		if(!std::isfinite(ph->energyTH)) {delete ph; continue;}
    }
}

//* 计算frameHessian的预计算值, 和状态的delta值
//@ 设置关键帧之间的关系
void FullSystem::setPrecalcValues()
{
    for(FrameHessian* fh : frameHessians)
    {
        fh->targetPrecalc.resize(frameHessians.size());
		for(unsigned int i=0;i<frameHessians.size();i++)  //? 还有自己和自己的???
			fh->targetPrecalc[i].set(fh, frameHessians[i], &Hcalib); // 计算Host 与 target之间的变换关系
    }

    ef->setDeltaF(&Hcalib);
}