#include "IOWrapper/ImageDisplay.h"
#include <iostream>

namespace IOWrap
{
std::unordered_set<std::string> openWindows;
boost::mutex openCVdisplayMutex;

void displayImage(const char* windowName, const cv::Mat& image, bool autoSize)
{
	if(disableAllDisplay) return;

	boost::unique_lock<boost::mutex> lock(openCVdisplayMutex);
	if(!autoSize)
	{
		if(openWindows.find(windowName) == openWindows.end())
		{
			cv::namedWindow(windowName, cv::WINDOW_NORMAL);
			cv::resizeWindow(windowName, image.cols, image.rows);
			openWindows.insert(windowName);
		}
	}
	cv::imshow(windowName, image);
}

void displayImage(const char* windowName, const MinimalImageB3* img, bool autoSize)
{
	displayImage(windowName, cv::Mat(img->h, img->w, CV_8UC3, img->data), autoSize);
}

}