#pragma once
#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/ximgproc.hpp"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace cv;
using namespace cv::ximgproc;

Mat binarizebw(Mat& src, double brightness, bool display_on);
Mat blobcleanup(Mat src_bin, int min_blob_area);
Mat histequalization(Mat bgr_image, bool display_on);
