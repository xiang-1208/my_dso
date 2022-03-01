#include "DatasetReader.h"


inline int getdir (std::string dir, std::vector<std::string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL)
    {
        return -1;
    }

    while ((dirp = readdir(dp)) != NULL) {
    	std::string name = std::string(dirp->d_name);

    	if(name != "." && name != "..")
    		files.push_back(name);
    }
    closedir(dp);


    std::sort(files.begin(), files.end());

    if(dir.at( dir.length() - 1 ) != '/') dir = dir+"/";
	for(unsigned int i=0;i<files.size();i++)
	{
		if(files[i].at(0) != '/')
			files[i] = dir + files[i];
		std::cout << files[i] << std::endl;
	}

    return files.size();
}

ImageFolderReader::ImageFolderReader(std::string path, std::string calibFile, std::string gammaFile, std::string vignetteFile)
: path(path),calibfile(calibFile),gammaFile(gammaFile),vignetteFile(vignetteFile)
{
    #if HAS_ZIPLIB
    	ziparchive=0;
		databuffer=0;
    #endif

    isZipped = (path.length()>4 && path.substr(path.length()-4) == ".zip");

	if(isZipped)
	{
#if HAS_ZIPLIB
		int ziperror=0;
		ziparchive = zip_open(path.c_str(),  ZIP_RDONLY, &ziperror);
		if(ziperror!=0)
		{
			printf("ERROR %d reading archive %s!\n", ziperror, path.c_str());
			exit(1);
		}
		files.clear();
		int numEntries = zip_get_num_entries(ziparchive, 0);
		for(int k=0;k<numEntries;k++)
		{
			const char* name = zip_get_name(ziparchive, k,  ZIP_FL_ENC_STRICT);
			std::string nstr = std::string(name);
			if(nstr == "." || nstr == "..") continue;
			files.push_back(name);
		}
		printf("got %d entries and %d files!\n", numEntries, (int)files.size());
		std::sort(files.begin(), files.end());
#else
		printf("ERROR: cannot read .zip archive, as compile without ziplib!\n");
		exit(1);
#endif
    }
	else
		getdir (path, files);


	undistort = Undistort::getUndistorterForFile(calibfile,gammaFile, vignetteFile);

	widthOrg = undistort->getOriginalSize()[0];
	heightOrg = undistort->getOriginalSize()[1];
	width=undistort->getSize()[0];
	height=undistort->getSize()[1];

	// load timestamps if possible.
	loadTimestamps();
	printf("ImageFolderReader: got %d files in %s!\n", (int)files.size(), path.c_str());
}

/********
 * 加载时间戳
 * 加载进私有变量，并且验证时间戳的有效性
*********/
void ImageFolderReader::loadTimestamps()
{
	std::ifstream tr;
	std::string timesFile = path.substr(0,path.find_last_of('/')) + "/times.txt";
	tr.open(timesFile.c_str());
	while(!tr.eof() && tr.good())
	{
		std::string line;
		char buf[1000];
		tr.getline(buf, 1000);

		int id;
		double stamp;
		float exposure = 0; //曝光时间

		if(sscanf(buf, "%d %lf %f", &id, &stamp, &exposure) == 3)
		{
			timestamps.push_back(stamp);
			exposures.push_back(exposure);
		}
		else if (sscanf(buf, "%d %lf", &id, &stamp) == 2)
		{
			timestamps.push_back(stamp);
			exposures.push_back(exposure);
		}
	}
	tr.close();

	// check if exposures are correct, (possibly skip)
	// 检查
	bool exposureisgood = ((int)exposures.size() == (int)getNumImage());
	for(int i=0;i<(int)exposures.size();i++)
	{
		if(exposures[i] == 0)
		{
			// 当曝光时间为0时进行修补
			float sum=0,num=0;
			if(i>0 && exposures[i-1] > 0) {sum += exposures[i-1]; num++;}
			if(i+1<(int)exposures.size() && exposures[i+1] > 0) {sum += exposures[i+1]; num++;}
			if(num>0)
				exposures[i] = sum/num;
		}
		if(exposures[i] == 0) exposureisgood=false;
	}

	if((int)getNumImage() != (int)timestamps.size())
	{
		printf("set timestamps and exposures to zero!\n");
		exposures.clear();
		timestamps.clear();
	}

	if((int)getNumImage() != (int)exposures.size() || !exposureisgood)
	{
		printf("set EXPOSURES to zero!\n");
		exposures.clear();
	}

	printf("got %d images and %d timestamps and %d exposures.!\n", (int)getNumImage(), (int)timestamps.size(), (int)exposures.size());
}

/*******
 * 
 * 
 * 
*******/
void ImageFolderReader::setGlobalCalibration()
{
	int w_out, h_out;
	Eigen::Matrix3f K;
	getCalibMono(K, w_out, h_out);
	setGlobalCalib(w_out, h_out, K);
}

void ImageFolderReader::getCalibMono(Eigen::Matrix3f &K, int &w, int &h)
{
	K = undistort->getK().cast<float>();
	// std::cout << K <<std::endl;
	w = undistort->getSize()[0];
	h = undistort->getSize()[1];
}

void ImageFolderReader::setGlobalCalib(int w,int h,Eigen::Matrix3f K)
{
	int wlvl = w;
	int hlvl = h;
	pyrLevelsUsed = 1;
	while(wlvl%2==0 && hlvl%2==0 && wlvl*hlvl > 5000 && pyrLevelsUsed < PYR_LEVELS)
	{
		wlvl /=2;
		hlvl /=2;
		pyrLevelsUsed++;
	}
	printf("using pyramid levels 0 to %d. coarsest resolution: %d x %d!\n",
			pyrLevelsUsed-1, wlvl, hlvl);
	if(wlvl>100 && hlvl > 100)
	{
		printf("\n\n===============WARNING!===================\n "
				"using not enough pyramid levels.\n"
				"Consider scaling to a resolution that is a multiple of a power of 2.\n");
	}
	if(pyrLevelsUsed < 3)
	{
		printf("\n\n===============WARNING!===================\n "
				"I need higher resolution.\n"
				"I will probably segfault.\n");
	}

	wM3G = w-3;
	hM3G = h-3;

	wG[0] = w;
	hG[0] = h;
	KG[0] = K;
	fxG[0] = K(0,0);
	fyG[0] = K(1,1);
	cxG[0] = K(0,2);
	cyG[0] = K(1,2);
	KiG[0] = KG[0].inverse();
	fxiG[0] = KiG[0](0,0);
	fyiG[0] = KiG[0](1,1);
	cxiG[0] = KiG[0](0,2);
	cyiG[0] = KiG[0](1,2);

	for (int level = 1; level < pyrLevelsUsed; ++ level)
	{
		wG[level] = w >> level;
		hG[level] = h >> level;

		fxG[level] = fxG[level-1] * 0.5;
		fyG[level] = fyG[level-1] * 0.5;

		cxG[level] = (cxG[0] + 0.5) / ((int)1<<level) - 0.5;
		cyG[level] = (cyG[0] + 0.5) / ((int)1<<level) - 0.5;

		KG[level]  << fxG[level], 0.0, cxG[level], 0.0, fyG[level], cyG[level], 0.0, 0.0, 1.0;	// synthetic
		KiG[level] = KG[level].inverse();
		
		fxiG[level] = KiG[level](0,0);
		fyiG[level] = KiG[level](1,1);
		cxiG[level] = KiG[level](0,2);
		cyiG[level] = KiG[level](1,2);
	}
}

double ImageFolderReader::getTimestamp(int id)
{
	if (timestamps.size() == 0) return id*0.1f;
	if (id > timestamps.size()) return 0;
	if (id < 0) return 0;
	return timestamps[id];
}

MinimalImageB* ImageFolderReader::getImageRaw_internal(int id, int unused)
{
	if (!isZipped)
	{
		return IOWrap::readImageBW_8U(files[id]);
	}
	else
	{
#if HAS_ZIPLIB
			if(databuffer==0) databuffer = new char[widthOrg*heightOrg*6+10000];
			zip_file_t* fle = zip_fopen(ziparchive, files[id].c_str(), 0);
			long readbytes = zip_fread(fle, databuffer, (long)widthOrg*heightOrg*6+10000);

			if(readbytes > (long)widthOrg*heightOrg*6)
			{
				printf("read %ld/%ld bytes for file %s. increase buffer!!\n", readbytes,(long)widthOrg*heightOrg*6+10000, files[id].c_str());
				delete[] databuffer;
				databuffer = new char[(long)widthOrg*heightOrg*30];
				fle = zip_fopen(ziparchive, files[id].c_str(), 0);
				readbytes = zip_fread(fle, databuffer, (long)widthOrg*heightOrg*30+10000);

				if(readbytes > (long)widthOrg*heightOrg*30)
				{
					printf("buffer still to small (read %ld/%ld). abort.\n", readbytes,(long)widthOrg*heightOrg*30+10000);
					exit(1);
				}
			}

			return (IOWrap::readStreamBW_8U(databuffer, readbytes));				
#else
			printf("ERROR: cannot read .zip archive, as compile without ziplib!\n");
			exit(1);
#endif	
	}
}

ImageAndExposure* ImageFolderReader::getImage_internal(int id, int unused)
{
	
	MinimalImageB* minimg = getImageRaw_internal(id, 0);
	ImageAndExposure* ret2 = undistort->undistort<unsigned char>(
			minimg,
			(exposures.size() == 0 ? 1.0f : exposures[id]),
			(timestamps.size() == 0 ? 0.0 : timestamps[id]));
	delete minimg;
	return ret2;
}