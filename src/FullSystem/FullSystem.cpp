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
            lock.unlock();
        }
    }
    else
    {
    }
}