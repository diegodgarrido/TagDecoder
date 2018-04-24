#include "globals.h"
#include "blobcleanup.h"

Mat binarizebw(Mat& src, double brightness, bool display_on)
{
    Mat mask,mask_gray,mask_er;
    int type = THRESH_BINARY;
    int blockSize = 29;
    double k_ = 8;
    double k = (double) ((k_ - 10.) / 10.);

    // be sure the mask is a gray image
    if (src.channels() == 3)
        cvtColor(src, mask_gray, COLOR_BGR2GRAY);
    else
        mask_gray = src.clone();

    mask_gray = 255 - mask_gray;
    if (display_on)
    {
        imshow("Inverse Gray", mask);
    }
    //adaptiveThreshold(mask,mask,255,CV_ADAPTIVE_THRESH_GAUSSIAN_C, CV_THRESH_BINARY,7,0);
    //threshold(mask, mask,25, 255, CV_THRESH_BINARY | CV_THRESH_OTSU);
    //threshold(mask, mask, 30, 255, CV_THRESH_BINARY );
    niBlackThreshold(mask_gray, mask, 255, type, blockSize, k);

    if ((brightness > BRIGHTNESS_LEVEL) && (brightness <= 1.0))
    {
        int erosion_elem = 0;
        int erosion_size = 2;
        int erosion_type = 0;

        if (erosion_elem == 0)
        {
            erosion_type = MORPH_RECT;
        }
        else if (erosion_elem == 1)
        {
            erosion_type = MORPH_CROSS;
        }
        else if (erosion_elem == 2)
        {
            erosion_type = MORPH_ELLIPSE;
        }


        Mat element = getStructuringElement(erosion_type,
                                            Size(2 * erosion_size + 1, 2 * erosion_size + 1),
                                            Point(erosion_size, erosion_size));
        /// Apply the erosion operation
        erode(mask, mask_er, element);
    }
    else
    {
        mask_er = mask;
    }

    return mask_er;
}



Mat blobcleanup(Mat src_bin, int min_blob_area)
{
    Mat src_clean_bin=src_bin;
    // find all contours in the binary image
    vector<vector<Point>> contours;
    vector<Vec4i> hierarchy;
    vector<int> small_blobs;
    findContours(src_clean_bin, contours, hierarchy, CV_RETR_CCOMP, CV_CHAIN_APPROX_SIMPLE);

    // Find indices of contours whose area is less than `threshold`
    if (!contours.empty())
    {
        for (size_t i = 0; i < contours.size(); ++i)
        {
            double contour_area = contourArea(contours[i]);
            if (contour_area < min_blob_area)
                small_blobs.push_back((int)i);
        }
    }

    // fill-in all small contours with zeros
    for (size_t i = 0; i < small_blobs.size(); ++i)
    {
        drawContours(src_clean_bin, contours, small_blobs[i], cv::Scalar(0), CV_FILLED, 8);
    }
    return(src_clean_bin);
}
