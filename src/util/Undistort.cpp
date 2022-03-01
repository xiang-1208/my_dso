#include "Undistort.h"


/**********
 * 
 * 
 * 
 * 
**************/
Undistort* Undistort::getUndistorterForFile(std::string configFilename, std::string gammaFilename, std::string vignetteFilename)
{
    printf ("Reading Calibration from file %s",configFilename.c_str());

    std::ifstream f(configFilename.c_str());
    if (!f.good())
    {
		f.close();
		printf(" ... not found. Cannot operate without calibration, shutting down.\n");
		f.close();
		return 0;
    }

    printf(" ... found!\n");
	std::string l1;
	std::getline(f,l1);
	f.close();

    float ic[10];

    Undistort* u;

    //RadTan (OpenCV) camera model
    if(std::sscanf(l1.c_str(), "%f %f %f %f %f %f %f %f",
			&ic[0], &ic[1], &ic[2], &ic[3],
			&ic[4], &ic[5], &ic[6], &ic[7]) == 8)
    {
        printf("found RadTan (OpenCV) camera model, building rectifier.\n");
        u = new UndistortRadTan(configFilename.c_str(), true);
		//if(!u->isValid()) {delete u; return 0; }        
    }
    //PINHOLE/ATAN camera model
    else if(std::sscanf(l1.c_str(), "%f %f %f %f %f",
			&ic[0], &ic[1], &ic[2], &ic[3], &ic[4]) == 5)
    {
		if(ic[4]==0)
		{
			printf("found PINHOLE camera model, building rectifier.\n");
            u = new UndistortPinhole(configFilename.c_str(), true);
			//if(!u->isValid()) {delete u; return 0; }
		}
		else
		{
			printf("found ATAN camera model, building rectifier.\n");
            u = new UndistortFOV(configFilename.c_str(), true);
			if(!u->isValid()) {delete u; return 0; }
		}    
    }
    else
    {
        printf("could not read calib file! exit.");
        exit(1);
    }
	u->loadPhotometricCalibration(gammaFilename,"",vignetteFilename);

    return u;
}




/**********
 * 载入光度校准
 * 
**********/
void Undistort::loadPhotometricCalibration(std::string file, std::string noiseImage, std::string vignetteImage)
{
	photometricUndist = new PhotometricUndistorter(file, noiseImage, vignetteImage,getOriginalSize()[0], getOriginalSize()[1]);
}

/***********
 * 
 * 
 * 
 * 
************/
void Undistort::readFromFile(const char* configFileName, int nPars, std::string prefix)
{
    photometricUndist=0;
	valid = false;
	passthrough=false;
	remapX = 0;
	remapY = 0;

    parsOrg = VecX(nPars);
    float outputCalibration[5];

	// read parameters
	std::ifstream infile(configFileName);
	assert(infile.good());

    std::string l1,l2,l3,l4;

	std::getline(infile,l1);
	std::getline(infile,l2);
    std::getline(infile,l3);
    std::getline(infile,l4);

    //l1 & l2
    if(nPars == 5) // fov model
    {
        char buf[1000];
        snprintf(buf, 1000, "%s%%lf %%lf %%lf %%lf %%lf", prefix.c_str());
        // std::cout<< "buf: " << buf << std::endl;
        // buf: %lf %lf %lf %lf %lf

        if(std::sscanf(l1.c_str(), buf, &parsOrg[0], &parsOrg[1], &parsOrg[2], &parsOrg[3], &parsOrg[4]) == 5 &&
				std::sscanf(l2.c_str(), "%d %d", &wOrg, &hOrg) == 2)
        {
			printf("Input resolution: %d %d\n",wOrg, hOrg);
			printf("In: %f %f %f %f %f\n",
					parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3], parsOrg[4]);            
        }
		else
		{
			printf("Failed to read camera calibration (invalid format?)\nCalibration file: %s\n", configFileName);
			infile.close();
			return;
		}
    }
    else if(nPars == 8) // KB, equi & radtan model
    {}
    else
    {}

    if(parsOrg[2] < 1 && parsOrg[3] < 1)
    {
        printf("\n\nFound fx=%f, fy=%f, cx=%f, cy=%f.\n I'm assuming this is the \"relative\" calibration file format,"
               "and will rescale this by image width / height to fx=%f, fy=%f, cx=%f, cy=%f.\n\n",
               parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3],
               parsOrg[0] * wOrg, parsOrg[1] * hOrg, parsOrg[2] * wOrg - 0.5, parsOrg[3] * hOrg - 0.5 );

        // rescale and substract 0.5 offset.
        // the 0.5 is because I'm assuming the calibration is given such that the pixel at (0,0)
        // contains the integral over intensity over [0,0]-[1,1], whereas I assume the pixel (0,0)
        // to contain a sample of the intensity ot [0,0], which is best approximated by the integral over
        // [-0.5,-0.5]-[0.5,0.5]. Thus, the shift by -0.5.
        parsOrg[0] = parsOrg[0] * wOrg;
        parsOrg[1] = parsOrg[1] * hOrg;
        parsOrg[2] = parsOrg[2] * wOrg - 0.5;
        parsOrg[3] = parsOrg[3] * hOrg - 0.5;
    }

    //l3
    std::cout << "l3: " << l3 << std::endl;
    std::cout << "outputCalibration[0]: " << outputCalibration[0] << std::endl;
    if(l3 == "crop")
    {}
    else if (l3 == "full")
    {}
    else if (std::sscanf(l3.c_str(), "%f %f %f %f %f", &outputCalibration[0], &outputCalibration[1], &outputCalibration[2], &outputCalibration[3], &outputCalibration[4]) == 5)
    {
		printf("Out: %f %f %f %f %f\n",
				outputCalibration[0], outputCalibration[1], outputCalibration[2], outputCalibration[3], outputCalibration[4]);
    }
    else
    {
		printf("Out: Failed to Read Output pars... not rectifying.\n");
		infile.close();
		return;        
    }

    //l4 
    if(std::sscanf(l4.c_str(),"%d %d", &w, &h)==2)
    {
        printf("Output resolution: %d %d\n",w, h);
    }
    else
    {
        printf("Out: Failed to Read Output resolution... not rectifying.\n");
		valid = false;
    }

	if(outputCalibration[0] == -1);
	else if(outputCalibration[0] == -2);
	else if(outputCalibration[0] == -3);
    else
    {
        if(outputCalibration[2] > 1 || outputCalibration[3] > 1)
        {
            printf("\n\n\nWARNING: given output calibration (%f %f %f %f) seems wrong. It needs to be relative to image width / height!\n\n\n",
                   outputCalibration[0],outputCalibration[1],outputCalibration[2],outputCalibration[3]);
        } 

        //按照输出大小重新设置内参
		K.setIdentity();
        K(0,0) = outputCalibration[0] * w;
        K(1,1) = outputCalibration[1] * h;
        K(0,2) = outputCalibration[2] * w - 0.5;
        K(1,2) = outputCalibration[3] * h - 0.5;    
    }

    remapX = new float[w*h];
    remapY = new float[w*h];

	for(int y=0;y<h;y++)
		for(int x=0;x<w;x++)
		{
			remapX[x+y*w] = x;
			remapY[x+y*w] = y;
		}

    distortCoordinates(remapX, remapY, remapX, remapY, h*w);

	for(int y=0;y<h;y++)
		for(int x=0;x<w;x++)
		{
			// make rounding resistant.
			float ix = remapX[x+y*w];
			float iy = remapY[x+y*w];

			if(ix == 0) ix = 0.001;
			if(iy == 0) iy = 0.001;
			if(ix == wOrg-1) ix = wOrg-1.001;
			if(iy == hOrg-1) ix = hOrg-1.001;

			if(ix > 0 && iy > 0 && ix < wOrg-1 &&  iy < wOrg-1)
			{
				remapX[x+y*w] = ix;
				remapY[x+y*w] = iy;
			}
			else
			{
				remapX[x+y*w] = -1;
				remapY[x+y*w] = -1;
			}
		}

	valid = true;

	printf("\nRectified Kamera Matrix:\n");
	std::cout << K << std::endl;
}

UndistortRadTan::UndistortRadTan(const char* configFileName, bool noprefix)
{
    printf("Creating RadTan undistorter\n");

    if(noprefix)
        readFromFile(configFileName, 8);
    else
        readFromFile(configFileName, 8,"RadTan");
}

UndistortPinhole::UndistortPinhole(const char* configFileName, bool noprefix)
{
    printf("Creating Pinhole undistorter\n");

    if(noprefix)
        readFromFile(configFileName, 5);
    else
        readFromFile(configFileName, 5,"Pinhole");
}

UndistortFOV::UndistortFOV(const char* configFileName, bool noprefix)
{
    printf("Creating FOV undistorter\n");

    if(noprefix)
        readFromFile(configFileName, 5);
    else
        readFromFile(configFileName, 5,"FOV");   
}


/****
 * FOV畸变模型
 * input:缩放图像每个坐标点
 * output：在原图上畸变后所在位置
******/
void UndistortFOV::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
//非静态成员函数后面加const（加到非成员函数或静态成员后面会产生编译错误），
//表示成员函数隐含传入的this指针为const指针，决定了在该成员函数中，
//任意修改它所在的类的成员的操作都是不允许的（因为隐含了对this指针的const引用）
{
    //std::cout << "distortCoordinates" << std::endl;
    float dist = parsOrg[4];
    float d2t = 2.0f*tan(dist/2.0f);

    //原始相机内参
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];

    //缩放相机内参
    float ofx = K(0,0);
	float ofy = K(1,1);
	float ocx = K(0,2);
	float ocy = K(1,2);

    //对缩放后的每一个像素
    for(int i=0;i<n;i++)
    {
        float x = in_x[i];  //缩放图像所在x坐标
        float y = in_y[i];  //缩放图像所在y坐标

        //归一化成像平面位置，原点在中央
        float ix = (x - ocx) / ofx;
		float iy = (y - ocy) / ofy;

        float r = sqrtf(ix*ix + iy*iy);
		float fac = (r==0 || dist==0) ? 1 : atanf(r * d2t)/(dist*r);

        //fac*ix 为畸变后的归一化平面x坐标
		ix = fx*fac*ix+cx;
		iy = fy*fac*iy+cy;

		out_x[i] = ix;
		out_y[i] = iy;
    }
}

void UndistortPinhole::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{}

void UndistortRadTan::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const
{}

PhotometricUndistorter::PhotometricUndistorter(
    std::string file, 
    std::string noiseImage, 
    std::string vignetteImage,
    int w_, int h_)
{
    valid = false;
	w = w_;
	h = h_;
    vignetteMap = 0;
	output = new ImageAndExposure(w,h);

	if(file=="" || vignetteImage=="")
	{
		printf("NO PHOTOMETRIC Calibration!\n");
	}

	// read G.
	std::ifstream f(file.c_str());
	printf("Reading Photometric Calibration from file %s\n",file.c_str());
	if (!f.good())
	{
		printf("PhotometricUndistorter: Could not open file!\n");
		return;
	}

    {
	    std::string line;
	    std::getline( f, line );
	    std::istringstream l1i( line );
	    std::vector<float> Gvec = std::vector<float>( std::istream_iterator<float>(l1i), std::istream_iterator<float>() ); 

        GDepth = Gvec.size();   //256

        if(GDepth < 256)
        {
            printf("PhotometricUndistorter: invalid format! got %d entries in first line, expected at least 256!\n",(int)Gvec.size());
            return;
        }

        for(int i=0;i<GDepth;i++) G[i] = Gvec[i];
        for(int i=0;i<GDepth-1;i++)
	    {
	    	if(G[i+1] <= G[i])
	    	{
	    		printf("PhotometricUndistorter: G invalid! it has to be strictly increasing, but it isnt!\n");
	    		return;
	    	}
	    }

	    float min=G[0];
        float max=G[GDepth-1];
        for(int i=0;i<GDepth;i++) G[i] = 255.0 * (G[i] - min) / (max-min);			// make it to 0..255 => 0..255.  
    }      


    printf("Reading Vignette Image from %s\n",vignetteImage.c_str());
    MinimalImage<unsigned short>* vm16 = IOWrap::readImageBW_16U(vignetteImage.c_str());
    MinimalImageB* vm8 = IOWrap::readImageBW_8U(vignetteImage.c_str());
    vignetteMap = new float[w*h];
    vignetteMapInv = new float[w*h];

    //dso:vm16
    if(vm16 != 0 )
    {
		if(vm16->w != w ||vm16->h != h)
		{
			printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
					vm16->w, vm16->h, w, h);
			if(vm16!=0) delete vm16;
			if(vm8!=0) delete vm8;
			return;
		}

        //最大值归一化
        float maxV=0;
        for (int i=0;i<w*h;i++)
            if(vm16->at(i)>maxV) maxV=vm16->at(i);
        for (int i=0;i<w*h;i++)
            vignetteMap[i] = vm16->at(i) / maxV;
    }
    else if (vm8 != 0 )
    {
		if(vm8->w != w ||vm8->h != h)
		{
			printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
					vm8->w, vm8->h, w, h);
			if(vm16!=0) delete vm16;
			if(vm8!=0) delete vm8;
			return;
		}

        //最大值归一化
        float maxV=0;
        for (int i=0;i<w*h;i++)
            if(vm8->at(i)>maxV) maxV=vm8->at(i);
        for (int i=0;i<w*h;i++)
            vignetteMap[i] = vm8->at(i) / maxV;       
    }
    else
    {
		printf("PhotometricUndistorter: Invalid vignette image\n");
		if(vm16!=0) delete vm16;
		if(vm8!=0) delete vm8;
		return;        
    }

	if(vm16!=0) delete vm16;
	if(vm8!=0) delete vm8;


	for(int i=0;i<w*h;i++)
		vignetteMapInv[i] = 1.0f / vignetteMap[i];

	printf("Successfully read photometric calibration!\n");
	valid = true;    
}


// 有问题
template<typename T>
ImageAndExposure* Undistort::undistort(const MinimalImage<T>* image_raw, float exposure, double timestamp, float factor) const
{
    if(image_raw->w != wOrg || image_raw->h != hOrg)
    {
		printf("Undistort::undistort: wrong image size (%d %d instead of %d %d) \n", image_raw->w, image_raw->h, w, h);
		exit(1);        
    }

    photometricUndist->processFrame<T>(image_raw->data, exposure, factor);
    ImageAndExposure* result = new ImageAndExposure(w, h, timestamp);
    photometricUndist->output->copyMetaTo(*result);

    if(!passthrough)
    {
        float* out_data = result->image;
        float* in_data = photometricUndist->output->image;

		float* noiseMapX=0;
		float* noiseMapY=0;

        if(benchmark_varNoise>0)
        {
            int numnoise=(benchmark_noiseGridsize+8)*(benchmark_noiseGridsize+8);
			noiseMapX=new float[numnoise];
			noiseMapY=new float[numnoise];
			memset(noiseMapX,0,sizeof(float)*numnoise);
			memset(noiseMapY,0,sizeof(float)*numnoise);

			for(int i=0;i<numnoise;i++)
			{
				noiseMapX[i] =  2*benchmark_varNoise * (rand()/(float)RAND_MAX - 0.5f);
				noiseMapY[i] =  2*benchmark_varNoise * (rand()/(float)RAND_MAX - 0.5f);
			}                        
        }

        for (int idx=w*h-1;idx>=0;idx--)
        {
            float xx = remapX[idx];
            float yy = remapY[idx];



			if(benchmark_varNoise>0)
			{
				float deltax = getInterpolatedElement11BiCub(noiseMapX, 4+(xx/(float)wOrg)*benchmark_noiseGridsize, 4+(yy/(float)hOrg)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );
				float deltay = getInterpolatedElement11BiCub(noiseMapY, 4+(xx/(float)wOrg)*benchmark_noiseGridsize, 4+(yy/(float)hOrg)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );
				float x = idx%w + deltax;
				float y = idx/w + deltay;
				if(x < 0.01) x = 0.01;
				if(y < 0.01) y = 0.01;
				if(x > w-1.01) x = w-1.01;
				if(y > h-1.01) y = h-1.01;

				xx = getInterpolatedElement(remapX, x, y, w);
				yy = getInterpolatedElement(remapY, x, y, w);
			}


            if(xx<0)
                out_data[idx] = 0;
            else
            {
                int xxi = xx;
                int yyi = yy;
                xx -= xxi;
                yy -= yyi;
                float xxyy = xx*yy;

                const float* src = in_data + xxi + yyi*wOrg;

                out_data[idx] = xxyy * src[1+wOrg]
									+ (yy-xxyy) * src[wOrg]
									+ (xx-xxyy) * src[1]
									+ (1-xx-yy+xxyy) * src[0];
            }
        }

		if(benchmark_varNoise>0)
		{
			delete[] noiseMapX;
			delete[] noiseMapY;
		}        
    }
	else
	{
		memcpy(result->image, photometricUndist->output->image, sizeof(float)*w*h);
	}  	

	applyBlurNoise(result->image);

	return result;    
}
template ImageAndExposure* Undistort::undistort<unsigned char>(const MinimalImage<unsigned char>* image_raw, float exposure, double timestamp, float factor) const;
template ImageAndExposure* Undistort::undistort<unsigned short>(const MinimalImage<unsigned short>* image_raw, float exposure, double timestamp, float factor) const;


void Undistort::applyBlurNoise(float* img) const
{
	if(benchmark_varBlurNoise==0) return;

	int numnoise=(benchmark_noiseGridsize+8)*(benchmark_noiseGridsize+8);
	float* noiseMapX=new float[numnoise];
	float* noiseMapY=new float[numnoise];
	float* blutTmp=new float[w*h];

	if(benchmark_varBlurNoise>0)
	{
		for(int i=0;i<numnoise;i++)
		{
				noiseMapX[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
				noiseMapY[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
		}
	}


	float gaussMap[1000];
	for(int i=0;i<1000;i++)
		gaussMap[i] = expf((float)(-i*i/(100.0*100.0)));

	// x-blur.
	for(int y=0;y<h;y++)
		for(int x=0;x<w;x++)
		{
			float xBlur = getInterpolatedElement11BiCub(noiseMapX,
					4+(x/(float)w)*benchmark_noiseGridsize,
					4+(y/(float)h)*benchmark_noiseGridsize,
					benchmark_noiseGridsize+8 );

			if(xBlur < 0.01) xBlur=0.01;


			int kernelSize = 1 + (int)(1.0f+xBlur*1.5);
			float sumW=0;
			float sumCW=0;
			for(int dx=0; dx <= kernelSize; dx++)
			{
				int gmid = 100.0f*dx/xBlur + 0.5f;
				if(gmid > 900 ) gmid = 900;
				float gw = gaussMap[gmid];

				if(x+dx>0 && x+dx<w)
				{
					sumW += gw;
					sumCW += gw * img[x+dx+y*this->w];
				}

				if(x-dx>0 && x-dx<w && dx!=0)
				{
					sumW += gw;
					sumCW += gw * img[x-dx+y*this->w];
				}
			}

			blutTmp[x+y*this->w] = sumCW / sumW;
		}

	// y-blur.
	for(int x=0;x<w;x++)
		for(int y=0;y<h;y++)
		{
			float yBlur = getInterpolatedElement11BiCub(noiseMapY,
					4+(x/(float)w)*benchmark_noiseGridsize,
					4+(y/(float)h)*benchmark_noiseGridsize,
					benchmark_noiseGridsize+8 );

			if(yBlur < 0.01) yBlur=0.01;

			int kernelSize = 1 + (int)(0.9f+yBlur*2.5);
			float sumW=0;
			float sumCW=0;
			for(int dy=0; dy <= kernelSize; dy++)
			{
				int gmid = 100.0f*dy/yBlur + 0.5f;
				if(gmid > 900 ) gmid = 900;
				float gw = gaussMap[gmid];

				if(y+dy>0 && y+dy<h)
				{
					sumW += gw;
					sumCW += gw * blutTmp[x+(y+dy)*this->w];
				}

				if(y-dy>0 && y-dy<h && dy!=0)
				{
					sumW += gw;
					sumCW += gw * blutTmp[x+(y-dy)*this->w];
				}
			}
			img[x+y*this->w] = sumCW / sumW;
		}

	delete[] noiseMapX;
	delete[] noiseMapY;
}

template<typename T>
void PhotometricUndistorter::processFrame(T* image_in, float exposure_time, float factor)
{
    int wh = w*h;
    float* data = output->image;
	assert(output->w == w && output->h == h);
	assert(data != 0);
	
	if(!valid || exposure_time <= 0 || setting_photometricCalibration==0)  
    {
        for(int i=0;i<wh;i++)
        {
            data[i] = factor*image_in[i];
        }
		output->exposure_time = exposure_time;
		output->timestamp = 0;
    } 
	else
	{
		for(int i=0; i<wh;i++)
		{
			data[i] = G[image_in[i]];
		}

		if(setting_photometricCalibration==2)
		{
			for(int i=0; i<wh;i++)
				data[i] *= vignetteMapInv[i];
		}

		output->exposure_time = exposure_time;
		output->timestamp = 0;
	}   
}
template void PhotometricUndistorter::processFrame<unsigned char>(unsigned char* image_in, float exposure_time, float factor);
template void PhotometricUndistorter::processFrame<unsigned short>(unsigned short* image_in, float exposure_time, float factor);

