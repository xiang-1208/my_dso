#pragma once
#include "util/MinimalImage.h"
#include <opencv2/highgui/highgui.hpp> 

namespace IOWrap
{
    MinimalImage<unsigned short>* readImageBW_16U(std::string filename);
    MinimalImageB* readImageBW_8U(std::string filename);
    MinimalImageB* readStreamBW_8U(char* data, int numBytes);
}