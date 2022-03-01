#pragma once
#include <string>
#include <fstream>
#include <iostream>

#include "Eigen/Core"
#include "NumType.h"
#include <iterator>
#include "IOWrapper/ImageRW.h"
#include "util/ImageAndExposure.h"
#include "util/globalutil.h"
#include "util/globalFuncs.h"

using namespace dso;

class PhotometricUndistorter
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    PhotometricUndistorter(std::string file, std::string noiseImage, std::string vignetteImage,int w_, int h_);
    template<typename T> void processFrame(T* image_in, float exposure_time, float factor=1);
    ~PhotometricUndistorter();

    ImageAndExposure* output;

    inline float* getG(){return G;};

private:
    bool valid; 
    int GDepth;
    float G[256*256];
    float* vignetteMap;
	float* vignetteMapInv;
    int w,h;

    
};

class Undistort
{
public:
    static Undistort* getUndistorterForFile(std::string configFilename, 
        std::string gammaFilename, 
        std::string vignetteFilename);

	template<typename T>
	ImageAndExposure* undistort(const MinimalImage<T>* image_raw, float exposure=0, double timestamp=0, float factor=1) const;
    
    void readFromFile(const char* configFileName, int nPars, std::string prefix = "");
    virtual void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const = 0;
    inline bool isValid() {return valid;};
    inline const Eigen::Vector2i getOriginalSize() {return Eigen::Vector2i(wOrg,hOrg);};
    inline const Eigen::Vector2i getSize() {return Eigen::Vector2i(w,h);};
    inline const Mat33 getK() const {return K;};
    void loadPhotometricCalibration(std::string file, std::string noiseImage, std::string vignetteImage);

    PhotometricUndistorter* photometricUndist;
    
protected:
    Mat33 K;
    int w, h,wOrg,hOrg;
    VecX parsOrg;
    bool valid;

    float* remapX;
    float* remapY;

    bool passthrough;

    void applyBlurNoise(float* img) const;
};

class UndistortFOV : public Undistort
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW; //eigen矩阵内存对齐

    UndistortFOV(const char* configFileName, bool noprefix);
    void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const;
};

class UndistortPinhole : public Undistort
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW; //eigen矩阵内存对齐

    UndistortPinhole(const char* configFileName, bool noprefix);
    void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const;
};

class UndistortRadTan : public Undistort
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW; //eigen矩阵内存对齐

    UndistortRadTan(const char* configFileName, bool noprefix);
    void distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const;
};