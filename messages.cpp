#include "opencv2/core.hpp"
#include "fileutils.h"
#include "globals.h"
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;
using namespace cv;

void TagMessage(const string message, string dirstring, string input_imageName, int index, ofstream &outputFile, string outtext, ifstream &goldFile, string goldtext, Tag rtag, string &output, float &confidence)
{
    confidence = 0;
    output = "empty";
    cout << "OCR_Tesseract  output = " << output << endl;
    cout << "Confidence :" << confidence << endl;
    cout << message << endl;
    string  id = std::to_string((long long)index);
    const string charsToRemove(".jpg");
    string file_name = input_imageName;
    removeStringFromString(file_name, charsToRemove, "");
    file_name += "_code.jpg";
    Mat empty_label = Mat::zeros(rtag.height_pxls, rtag.width_pxls, CV_8UC1);
    imwrite(file_name, empty_label);
    if (!outtext.empty())
    {
        if (!dirstring.empty())
        {
            string short_file_name = file_name;

            removeStringFromString(short_file_name, dirstring, "");
            if (!goldtext.empty())
            {
                string  gfile_name, label_code;
                goldFile >> gfile_name >> label_code;
                outputFile << setw(40) << short_file_name << "   " << setw(6) <<  output << "   " << confidence << "   " << "Gold: " << setw(6) << label_code << endl;
                //outputFile << print("%40s    %6s    %2.2f    Gold: %6s\n", short_file_name.c_str(), output.c_str(), confidence, label_code.c_str());
            }
            else
            {
                outputFile << setw(40) << short_file_name << "   " << setw(6) <<  output << "   " << confidence << endl;
                //outputFile << print("%40s    %6s    %2.2f\n", short_file_name.c_str(), output.c_str(), confidence);
            }
        }
        else
        {
            outputFile << setw(40) << file_name << "   " << setw(6) <<  output << "   " << confidence << endl;
            //outputFile << print("%40s    %6s    %2.2f\n", file_name.c_str(), output.c_str(), confidence);
            //outputFile << file_name << " " << output << " " << confidence << endl;
        }
    }
    return;
}
void WriteMessage(const Mat final_mask,int &correct_labels, string dirstring, string input_imageName, int index, ofstream &outputFile, string outtext, ifstream &goldFile, string goldtext, Tag rtag, string &output, float &confidence)
{
    string  id = std::to_string((long long)index);
    const string charsToRemove(".jpg");
    string file_name = input_imageName;

    removeStringFromString(file_name, charsToRemove, "");
    file_name += "_code.jpg";
    imwrite(file_name, final_mask);
    if (!outtext.empty())
    {
        if (!dirstring.empty())
        {
            string short_file_name = file_name;
            removeStringFromString(short_file_name, dirstring, "");
            if (!goldtext.empty())
            {
                string  gfile_name, label_code;
                goldFile >> gfile_name >> label_code;
                outputFile << setw(40) << short_file_name << "   " << setw(6) <<  output << "   " << confidence << "   " << "Gold: " << setw(6) << label_code << endl;
                //outputFile << print("%40s    %6s    %2.2f    Gold: %6s\n", short_file_name.c_str(), output.c_str(), confidence, label_code.c_str());
                if (label_code == output)
                {
                    correct_labels++;
                }
            }
            else
            {
                outputFile << setw(40) << short_file_name << "   " << setw(6) <<  output << "   " << confidence << endl;
                //outputFile << print("%40s    %6s    %2.2f\n", short_file_name.c_str(), output.c_str(), confidence);
            }
        }
        else
        {
            outputFile << setw(40) << file_name << "   " << setw(6) <<  output << "   " << confidence << endl;
            //outputFile << print("%40s    %6s    %2.2f\n", file_name.c_str(), output.c_str(), confidence);
        }
    }
}
