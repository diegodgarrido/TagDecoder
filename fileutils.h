#pragma once
#ifndef __GNUC__
#include <Windows.h>
#endif

#include "opencv2/core.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef __GNUC__
#include <dirent.h>
#else
#include "dirent.h"
#endif


#include <string>
#include <iostream>
#include <vector>




using namespace std;
//void lsfiles(string folder, vector<string> &files);
int getFrameSuffixes(const char *dirname, const string suffix_name, vector<string> &suffix_names);
int getFramePrefixes(const char *dirname, const string prefix_name, vector<string> &prefix_names);
bool getFrameNames(const string prefix_name, const vector<string> &prefix_names, string index_name_begin, string index_name_end, vector<string> &frame_names);
bool removeStringFromString(std::string& str, const std::string& from, const std::string& to);
