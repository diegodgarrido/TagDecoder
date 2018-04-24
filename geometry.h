#pragma once
#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include "globals.h"


using namespace std;
using namespace cv;


int CircularDetection(Mat &img,bool contrast_stretch, const int min_dist_center, const int min_number_of_votes, const int min_radius, const int max_radius, vector<Vec3f> &circles);
int CircularDetection(Mat &img,const int min_dist_center, const int min_number_of_votes, const int min_radius, const int max_radius, vector<Vec3f> &circles);
int ConvertQuadinRects(const Mat src, const vector<Vec3f>& Quadrilateral, Rect & Rectangle);
int ConvertQuadinLabelBounds(const Mat src, const vector<Vec3f>& Quadrilateral, Rect & Rectangle);
bool FindQuad4pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad, int &scale);
bool FindQuad3pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad, int &scale);
bool isValidQuad(vector<Vec3f> &Quad, double qwidth, double qheight, int tol);
Rect IncreaseRect(Mat src, Rect R, double val);
Rect ReduceRect(Rect R, double val);
double CalcAngle(Mat src_gray, const vector<Vec3f>& Quadrilateral, double min_angle_to_assign_zero, double max_angle_to_assign_zero);
Mat Rotate(Mat src, double angle);
void PrintCircles(vector<Vec3f> &circles);
void DrawCircles(const Mat image, vector<Vec3f> &circles, string name, double scale);
void ShowFrame(const Mat image, string name, double scale);


