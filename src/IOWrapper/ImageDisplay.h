#pragma once
#include "util/MinimalImage.h"
#include <boost/thread.hpp>
#include "util/globalutil.h"
#include <unordered_set>
#include <opencv2/highgui/highgui.hpp> 

using namespace dso;

namespace IOWrap
{
void displayImage(const char* windowName, const MinimalImageB3* img, bool autoSize = false);
}