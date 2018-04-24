#pragma once
#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/photo/photo.hpp"

#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace cv;

Mat histequalization(Mat src, bool display_on);
Mat equalizeIntensity(const Mat& inputImage);
Mat ConstrastStretching(const Mat src, int r1, int s1, int r2, int s2);
Mat& GammaCorrection(Mat& I, double fGamma);
void getImgPerception(const Mat& src, double &brightness, double &contrast, double &saturation);
Mat DecreaseLuminance(Mat& src, int val);
double EstimateAngle(Mat src_gray, int threshold, double minLineLenght, double maxLineGap, double min_angle_to_assign_zero, double max_angle_to_assign_zero, bool display_on);
