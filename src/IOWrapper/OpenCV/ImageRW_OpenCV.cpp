#include "IOWrapper/ImageRW.h"
#include <iostream>

namespace IOWrap
{
MinimalImage<unsigned short>* readImageBW_16U(std::string filename)
{
	cv::Mat m = cv::imread(filename, CV_LOAD_IMAGE_UNCHANGED);
	if(m.rows*m.cols==0)
	{
		printf("cv::imread could not read image %s! this may segfault. \n", filename.c_str());
		return 0;
	}
	if(m.type() != CV_16U)
	{
		printf("readImageBW_16U called on image that is not a 16bit grayscale image. this may segfault. \n");
		return 0;
	}
	MinimalImage<unsigned short>* img = new MinimalImage<unsigned short>(m.cols, m.rows);
	memcpy(img->data, m.data, 2*m.rows*m.cols);
	return img;
}

MinimalImageB* readImageBW_8U(std::string filename)
{
	cv::Mat m = cv::imread(filename, CV_LOAD_IMAGE_UNCHANGED);
	if(m.rows*m.cols==0)
	{
		printf("cv::imread could not read image %s! this may segfault. \n", filename.c_str());
		return 0;
	}
	if(m.type() != CV_8U)
	{
		printf("readImageBW_8U called on image that is not a 8bit grayscale image. this may segfault. \n");
		return 0;
	}
	MinimalImageB* img = new MinimalImageB(m.cols, m.rows);
	memcpy(img->data, m.data, 2*m.rows*m.cols);
	return img;
}

MinimalImageB* readStreamBW_8U(char* data, int numBytes)
{
	cv::Mat m = cv::imdecode(cv::Mat(numBytes,1,CV_8U, data), CV_LOAD_IMAGE_GRAYSCALE);
	if(m.rows*m.cols==0)
	{
		printf("cv::imdecode could not read stream (%d bytes)! this may segfault. \n", numBytes);
		return 0;
	}
	if(m.type() != CV_8U)
	{
		printf("readImageBW_8U called on image that is not a 8bit grayscale image. this may segfault. \n");
		return 0;
	}
	MinimalImageB* img = new MinimalImageB(m.cols, m.rows);
	memcpy(img->data, m.data, m.rows*m.cols);
	return img;
}

}