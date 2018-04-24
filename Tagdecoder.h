#ifndef TAGDECODER_H_INCLUDED
#define TAGDECODER_H_INCLUDED
#ifdef __GNUC__
#include "tesseract/baseapi.h" // Tesseract header
#else
#include "baseapi.h"
#endif
#include "fileutils.h"
#include "geometry.h"
#include "blobcleanup.h"
#include "imageproc.h"
#include <iostream>
#include <fstream>
#include <string>

void TagMessage(const string message, string dirstring, string input_imageName, int index,  ofstream &outputFile, string outtext, ifstream &goldFile, string goldtext, Tag rtag, string &output, float &confidence);
void WriteMessage(const Mat final_mask, int &correct_labels, string dirstring, string input_imageName, int index, ofstream &outputFile, string outtext, ifstream &goldFile, string goldtext, Tag rtag, string &output, float &confidence);
using namespace std;
#endif // TAGDECODER_H_INCLUDED
