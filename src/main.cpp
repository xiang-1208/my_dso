#include <stdio.h>


// std::cout << ": " <<  << std::endl;
#include <thread>
#include "util/DatasetReader.h"
#include "util/ImageAndExposure.h"
#include "FullSystem/FullSystem.h"
#include "IOWrapper/Output3DWrapper.h"
#include "IOWrapper/Pangolin/PangolinDSOViewer.h"



void parseArgument(char* arg);

int main(int argc, char** argv )
{
    //参数读取
    for(int i=1; i<argc;i++)
	    parseArgument(argv[i]);

    //畸变文件读取
	//xiang finish 2022.1.20
    ImageFolderReader* reader = new ImageFolderReader(source,calib, gammaCalib, vignette);
	//各层金字塔外参预计算
	//xiang finish 2022.1.24
	reader->setGlobalCalibration();

	if(setting_photometricCalibration>0 && reader->getPhotometricGamma() == 0 )
	{
		printf("ERROR: dont't have photometric calibation. Need to use commandline options mode=1 or mode=2 ");
		exit(1);		
	}

	int lstart=start;
	int lend = end;
	int linc = 1;
	if (reverse)
	{
		lstart = end-1;
		if (lstart > reader->getNumImages())
			lstart = reader->getNumImages()-1;
		lend = start;
		linc = -1;
	}

	FullSystem* fullSystem = new FullSystem();
	//???
	//xiang finish 2022.1.25
	fullSystem->setGammaFunction(reader->getPhotometricGamma());
	fullSystem->linearizeOperation = (playbackSpeed==0);

	//显示模块初始化
	//xiang finish 2022.2.6
    IOWrap::PangolinDSOViewer* viewer = 0;
	if(!disableAllDisplay)
    {
        viewer = new IOWrap::PangolinDSOViewer(wG[0],hG[0], false);
        fullSystem->outputWrapper.push_back(viewer);
    }

	std::thread runthread([&]() {
        std::vector<int> idsToPlay;
        std::vector<double> timesToPlayAt;	
		for(int i=lstart;i>= 0 && i< reader->getNumImages() && linc*i < linc*lend;i+=linc)	
		{
			idsToPlay.push_back(i);
			if (timesToPlayAt.size()==0)
				timesToPlayAt.push_back((double)0);
			else
			{
				double tsThis = reader->getTimestamp(idsToPlay[idsToPlay.size()-1]);
				double tsPrev = reader->getTimestamp(idsToPlay[idsToPlay.size()-2]);
				timesToPlayAt.push_back(timesToPlayAt.back() +  fabs(tsThis-tsPrev)/playbackSpeed);
			}
		}
		
		// std::vector<ImageA

		struct timeval tv_start;
		gettimeofday(&tv_start,NULL);
        clock_t started = clock();
        double sInitializerOffset=0;

		for(int ii=0;ii<(int)idsToPlay.size(); ii++)
		{
			if(!fullSystem->initialized)
			{
				gettimeofday(&tv_start, NULL);
                started = clock();
                sInitializerOffset = timesToPlayAt[ii];
			}
			int i = idsToPlay[ii];

			ImageAndExposure* img;
            if(preload)
			{
                //img = preloadedImages[ii];
			}
            else
			{
				//读取图像并且通过校准
				//xiang finish 2022.2.14
                img = reader->getImage(i);
			}

			bool skipFrame=false;

			if(playbackSpeed!=0)
			{
				struct timeval tv_now; 
				gettimeofday(&tv_now, NULL);
				double sSinceStart = sInitializerOffset + ((tv_now.tv_sec-tv_start.tv_sec) + (tv_now.tv_usec-tv_start.tv_usec)/(1000.0f*1000.0f));
 
                if(sSinceStart < timesToPlayAt[ii])
                    usleep((int)((timesToPlayAt[ii]-sSinceStart)*1000*1000));
                else if(sSinceStart > timesToPlayAt[ii]+0.5+0.1*(ii%2))
                {
                    printf("SKIPFRAME %d (play at %f, now it is %f)!\n", ii, timesToPlayAt[ii], sSinceStart);
                    skipFrame=true;
                }
			}

			std::cout << "image of " << i <<std::endl;
			if(!skipFrame) fullSystem->addActiveFrame(img, i);

			delete img;

			if(fullSystem->initFailed || setting_fullResetRequested)
			{
				if(ii < 250 || setting_fullResetRequested)
                {
                	printf("RESETTING!\n");

                	std::vector<IOWrap::Output3DWrapper*> wraps = fullSystem->outputWrapper;
                	delete fullSystem;

                	for(IOWrap::Output3DWrapper* ow : wraps) ow->reset();

                	fullSystem = new FullSystem();
                	fullSystem->setGammaFunction(reader->getPhotometricGamma());
                	fullSystem->linearizeOperation = (playbackSpeed==0);


                	fullSystem->outputWrapper = wraps;

                	setting_fullResetRequested=false;	
				}			
			}

            if(fullSystem->isLost)
            {
                    printf("LOST!!\n");
                    break;
            }
		}

        clock_t ended = clock();
        struct timeval tv_end;
        gettimeofday(&tv_end, NULL);	

        int numFramesProcessed = abs(idsToPlay[0]-idsToPlay.back());
        double numSecondsProcessed = fabs(reader->getTimestamp(idsToPlay[0])-reader->getTimestamp(idsToPlay.back()));
        double MilliSecondsTakenSingle = 1000.0f*(ended-started)/(float)(CLOCKS_PER_SEC);
        double MilliSecondsTakenMT = sInitializerOffset + ((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f);
        printf("\n======================"
                "\n%d Frames (%.1f fps)"
                "\n%.2fms per frame (single core); "
                "\n%.2fms per frame (multi core); "
                "\n%.3fx (single core); "
                "\n%.3fx (multi core); "
                "\n======================\n\n",
                numFramesProcessed, numFramesProcessed/numSecondsProcessed,
                MilliSecondsTakenSingle/numFramesProcessed,
                MilliSecondsTakenMT / (float)numFramesProcessed,
                1000 / (MilliSecondsTakenSingle/numSecondsProcessed),
                1000 / (MilliSecondsTakenMT / numSecondsProcessed));

        if(setting_logStuff)
        {
            std::ofstream tmlog;
            tmlog.open("logs/time.txt", std::ios::trunc | std::ios::out);
            tmlog << 1000.0f*(ended-started)/(float)(CLOCKS_PER_SEC*reader->getNumImages()) << " "
                  << ((tv_end.tv_sec-tv_start.tv_sec)*1000.0f + (tv_end.tv_usec-tv_start.tv_usec)/1000.0f) / (float)reader->getNumImages() << "\n";
            tmlog.flush();
            tmlog.close();
        }				
	});

    if(viewer != 0)
        viewer->run();

	runthread.join();

	for(IOWrap::Output3DWrapper* ow : fullSystem->outputWrapper)
	{
		ow->join();
		delete ow;
	}

	printf("DELETE FULLSYSTEM!\n");
	delete fullSystem;

	printf("DELETE READER!\n");
	delete reader;

	printf("EXIT NOW!\n");
	return 0;
	return 0;
}




void parseArgument(char* arg)
{
    int option;
    float foption;
	char buf[1000];

	if(1==sscanf(arg,"start=%d",&option))
	{
		start = option;
		printf("START AT %d!\n",start);
		return;
	}

	if(1==sscanf(arg,"end=%d",&option))
	{
		end = option;
		printf("END AT %d!\n",start);
		return;
	}

	if(1==sscanf(arg,"files=%s",buf))
	{
		source = buf;
		printf("loading data from %s!\n", source.c_str());
		return;
	}

	if(1==sscanf(arg,"calib=%s",buf))
	{
		calib = buf;
		printf("loading calibration from %s!\n", calib.c_str());
		return;
	}

	if(1==sscanf(arg,"undistort=%s",buf))
	{
		undistortion = buf;
		printf("loading calibration from %s!\n", undistortion.c_str());
		return;
	}

	if(1==sscanf(arg,"extrinsic=%s",buf))
	{
		extrinsic = buf;
		printf("loading extrinsic from %s!\n", extrinsic.c_str());
		return;
	}

	if(1==sscanf(arg,"vignette=%s",buf))
	{
		vignette = buf;
		printf("loading vignette from %s!\n", vignette.c_str());
		return;
	}

	if(1==sscanf(arg,"gamma=%s",buf))
	{
		gammaCalib = buf;
		printf("loading gammaCalib from %s!\n", gammaCalib.c_str());
		return;
	}

	if(1==sscanf(arg,"reverse=%d",&option))
	{
		if (option = 1)
		{
			reverse = true;
			printf("REVERSE!\n");
		}
		return;
	}

	if(1==sscanf(arg,"nolog=%d",&option))
	{
		if(option==1)
		{
			setting_logStuff = false;
			printf("DISABLE LOGGING!\n");
		}
		return;
	}

	if(1==sscanf(arg,"speed=%f",&foption))
	{
		playbackSpeed = foption;
		printf("PLAYBACK SPEED %f!\n", playbackSpeed);
		return;
	}

	if(1==sscanf(arg,"nogui=%d",&option))
	{
		if(option==1)
		{
			disableAllDisplay = true;
			printf("NO GUI!\n");
		}
		return;
	}

	if(1==sscanf(arg,"mode=%d",&option))
	{
		mode = option;
		if(option==0)
		{
			printf("PHOTOMETRIC MODE WITH CALIBRATION!\n");
		}
		if(option==1)
		{
			printf("PHOTOMETRIC MODE WITHOUT CALIBRATION!\n");
			setting_photometricCalibration = 0;
			setting_affineOptModeA = 0; //-1: fix. >=0: optimize (with prior, if > 0).
			setting_affineOptModeB = 0; //-1: fix. >=0: optimize (with prior, if > 0).
		}
		if(option==2)
		{
			printf("PHOTOMETRIC MODE WITH PERFECT IMAGES!\n");
			setting_photometricCalibration = 0;
			setting_affineOptModeA = -1; //-1: fix. >=0: optimize (with prior, if > 0).
			setting_affineOptModeB = -1; //-1: fix. >=0: optimize (with prior, if > 0).
            setting_minGradHistAdd=3;
		}
		return;
	}

	printf("could not parse argument \"%s\"!!!!\n", arg);
}
