#pragma once
#include <string>
#if HAS_ZIPLIB
	#include "zip.h"
    #include "zlib.h"
#endif
#include <iostream>
#include <vector>
#include <algorithm>
#include <dirent.h>
#include "Undistort.h"
#include <Eigen/Dense>
#include "util/globalutil.h"
#include "util/ImageAndExposure.h"
#include "IOWrapper/ImageRW.h"

using namespace dso;

class ImageFolderReader
{
public:
    ImageFolderReader(std::string path, std::string calibFile, std::string gammaFile, std::string vignetteFile);
    void setGlobalCalibration();

    double getTimestamp(int id);

    inline float* getPhotometricGamma()
	{
		if(undistort==0 || undistort->photometricUndist==0) return 0;
		return undistort->photometricUndist->getG();
	};

    inline ImageAndExposure* getImage(int id, bool forceLoadDirectly=false)
	{
		return getImage_internal(id, 0);
	}

    inline int getNumImages()
	{
		return files.size();
	};

private:
    std::string path;
    std::string calibfile;
    std::string gammaFile;
    std::string vignetteFile;
    std::vector<std::string> files;

    int width, height;
	int widthOrg, heightOrg;

    bool isZipped;
    zip_t* ziparchive;
	char* databuffer;
    Undistort* undistort;
    std::vector<double> timestamps;
    std::vector<float> exposures;

    MinimalImageB* getImageRaw_internal(int id, int unused);
    ImageAndExposure* getImage_internal(int id, int unused);
    void getCalibMono(Eigen::Matrix3f &K, int &w, int &h);
    void loadTimestamps();
    void setGlobalCalib(int w,int h,Eigen::Matrix3f K);
    inline int getNumImage() {return (files.size());};
};
