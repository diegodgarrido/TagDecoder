
#include "imageproc.h"
Mat histequalization(Mat bgr_image,bool display_on)
{
    // READ RGB color image and convert it to Lab
    Mat lab_image;
    cvtColor(bgr_image, lab_image, CV_BGR2Lab);

    // Extract the L channel
    std::vector<Mat> lab_planes(3);
    split(lab_image, lab_planes);  // now we have the L image in lab_planes[0]

    // apply the CLAHE algorithm to the L channel
    Ptr<CLAHE> clahe = createCLAHE();
    clahe->setClipLimit(4);
    Mat dst;
    clahe->apply(lab_planes[0], dst);

    // Merge the the color planes back into an Lab image
    dst.copyTo(lab_planes[0]);
    merge(lab_planes, lab_image);

    // convert back to RGB
    Mat image_clahe;
    cvtColor(lab_image, image_clahe, CV_Lab2BGR);

    // display the results  (you might also want to see lab_planes[0] before and after).
    if (display_on)
    {
        Mat resized;
        resize(bgr_image, resized, Size(bgr_image.cols >> 2, bgr_image.rows >> 2));
        imshow("image original", resized);
        resize(image_clahe, resized, Size(bgr_image.cols >> 2, bgr_image.rows >> 2));
        imshow("image CLAHE", resized);
    }
    return(image_clahe);
}
Mat equalizeIntensity(const Mat& inputImage)
{
    if (inputImage.channels() >= 3)
    {
        Mat ycrcb;

        cvtColor(inputImage, ycrcb, CV_BGR2YCrCb);

        vector<Mat> channels;
        split(ycrcb, channels);

        equalizeHist(channels[0], channels[0]);

        Mat result;
        merge(channels, ycrcb);

        cvtColor(ycrcb, result, CV_YCrCb2BGR);

        return result;
    }
    return Mat();
}

Mat& GammaCorrection(Mat& I, double fGamma)
{
    CV_Assert(I.data);

    // accept only char type matrices
    CV_Assert(I.depth() != sizeof(uchar));

    // build look up table
    unsigned char lut[256];
    for (int i = 0; i < 256; i++)
    {
        lut[i] = (unsigned char)(pow((double)(i / 255.0), 1. / fGamma) * 255.0);
    }

    const int channels = I.channels();
    switch (channels)
    {
    case 1:
    {

        MatIterator_<uchar> it, end;
        for (it = I.begin<uchar>(), end = I.end<uchar>(); it != end; it++)
            //*it = pow((float)(((*it))/255.0), fGamma) * 255.0;
            *it = lut[(*it)];

        break;
    }
    case 3:
    {

        MatIterator_<Vec3b> it, end;
        for (it = I.begin<Vec3b>(), end = I.end<Vec3b>(); it != end; it++)
        {
            //(*it)[0] = pow((float)(((*it)[0])/255.0), fGamma) * 255.0;
            //(*it)[1] = pow((float)(((*it)[1])/255.0), fGamma) * 255.0;
            //(*it)[2] = pow((float)(((*it)[2])/255.0), fGamma) * 255.0;
            (*it)[0] = lut[((*it)[0])];
            (*it)[1] = lut[((*it)[1])];
            (*it)[2] = lut[((*it)[2])];
        }

        break;
    }
    }

    return I;
}

int computeOutput(int x, int r1, int s1, int r2, int s2)
{
    float result=((float) x);
    float fx = (float)x;
    float fr1 = (float)r1;
    float fs1 = (float)s1;
    float fr2 = (float)r2;
    float fs2 = (float)s2;

    if (0 <= x && x <= r1)
    {
        result = round(fs1 / fr1 * fx);
    }
    else if (r1 < x && x <= r2)
    {
        result = round(((fs2 - fs1) / (fr2 - fr1)) * (fx - fr1) + fs1);
    }
    else if (r2 < x && x <= 255)
    {
        result = round(((255.f - fs2) / (255.f - fr2)) * (fx - fr2) + fs2);
    }
    return ((int)result);
}

Mat ConstrastStretching(const Mat src, int r1, int s1, int r2, int s2)
{
    Mat img = src.clone();
    int cc;
    if (src.channels() == 3)
        cc = 3;
    else
        cc = 1;

    if (cc == 3)
    {
        for (int y = 0; y < img.rows; y++)
        {
            for (int x = 0; x < img.cols; x++)
            {
                for (int c = 0; c < cc; c++)
                {
                    int output = computeOutput(src.at<Vec3b>(y, x)[c], r1, s1, r2, s2);
                    img.at<Vec3b>(y, x)[c] = saturate_cast<uchar>(output);

                }
            }
        }
    }
    else
    {
        for (int y = 0; y < img.rows; y++)
        {
            for (int x = 0; x < img.cols; x++)
            {
                int output = computeOutput(src.at<uchar>(y, x), r1, s1, r2, s2);
                img.at<uchar>(y, x) = saturate_cast<uchar>(output);
            }
        }
    }
    return(img);
}

void getImgPerception(const Mat& src, double &brightness, double &contrast, double &saturation)
{
    Mat img;
    double minVal;
    double maxVal;
    Point minLoc;
    Point maxLoc;

    if (src.channels() == 3)
        cvtColor(src, img, CV_BGR2GRAY);
    else
        img = src.clone();

    /// Sum all luminance pixels
    Scalar summ = sum(img);
    //-- percentage conversion factor brightness [0-1.0]
    brightness = summ[0] / ((pow(2, 8) - 1)*img.rows * img.cols);
#if 0
    /// Establish the number of bins
    float range[] = { 0, 256 };
    int histSize = 256;
    const float* histRange = { range };
    bool uniform = true;
    bool accumulate = false;
    Mat hist;

    /// Compute the histograms:

    calcHist(&img, 1, 0, Mat(), hist, 1, &histSize, &histRange, uniform, accumulate);
    hist /= (double) img.total();

    double entropy = 0.;
    for (int i = 0; i < 256; i++)
    {
        float pi = hist.at<float>(i, 0);
        if (pi > 0)
        {
            entropy += -pi * log2(pi);
        }
    }
    //-- percentage conversion factor brightness [0-8.0]
    contrast = entropy; // 1/8 *
#endif
    // Calculate RMS contrast

    contrast = 0.;

    for (int y = 0; y < src.rows; y++)
    {
        for (int x = 0; x < src.cols; x ++)
        {

            contrast += pow(((src.at<unsigned char>(y, x) / 255.0) - brightness), 2.0);
        }
    }

    contrast /= ((float)src.rows) * ((float)src.rows);
    contrast = pow(contrast, 0.5f);

    //-- saturation
    minMaxLoc(img, &minVal, &maxVal, &minLoc, &maxLoc);
    saturation = (maxVal - minVal) / (maxVal + minVal);


    return;
}


Mat DecreaseLuminance(Mat& src, int val)
{

    Mat src_hsv, src_bgr;

    // BGR to HSV
    cvtColor(src, src_hsv, CV_BGR2HSV);

    for (int i=0; i < src.rows ; i++)
    {
        for(int j=0; j < src.cols; j++)
        {
            // You need to check this, but I think index 1 is for saturation, but it might be 0 or 2
            int idx = 1;
            int pxlv = (int) src_hsv.at<cv::Vec3b>(i,j)[idx];
            //cout << pxlv << endl;
            pxlv = pxlv - val;
            int a = ((pxlv < 0) ? 0: pxlv);
            src_hsv.at<cv::Vec3b>(i,j)[idx] = (uchar)a;
            // or:
            // img.at<cv::Vec3b>(i,j)[idx] += adds_constant_value;
        }
    }

    // HSV back to BGR
    cvtColor(src_hsv, src_bgr, CV_HSV2BGR);
    return(src_bgr);
}

static double makehist(vector<double> samples, int histsize, double maxval, double minval)
{
    int i = 0, j = 0;
    vector<int> vFreq(histsize);
    vector<double> vBins(histsize);
    double histBinLow, histBinHigh, Sample;
    double delta = (maxval - minval) / ((double)histsize);

    for (int k = 0; k < histsize; k++)
    {
        vFreq[k] = 0;
        vBins[k] = 0;
    }
    i = 0;
    while (i < (int)samples.size())
    {
        histBinLow = minval;
        j = 0;
        while (j < histsize)
        {
            histBinHigh = histBinLow + delta;
            Sample = samples.at(i);
            if (Sample >= histBinLow && Sample < histBinHigh)
            {
                vFreq[j] = vFreq[j] + 1;
                vBins[j] += Sample;
                break;
            }
            j++;
            histBinLow = histBinHigh;
        }
        i++;
    }

    int maxFreq = 0;
    int maxidx = 0;
    for (int k = 0; k < (int)vFreq.size(); k++)
    {

        if (vFreq[k] > maxFreq)
        {
            maxFreq = vFreq[k];
            maxidx = k;
        }
    }
    double angle = vBins[maxidx] / ((double)vFreq[maxidx]);

    return(angle);

}
double EstimateAngle(Mat src_gray, int threshold, double minLineLenght, double maxLineGap, double min_angle_to_assign_zero, double max_angle_to_assign_zero, bool display_on)
{

    Mat src0 = src_gray;
    bitwise_not(src0, src0);
    vector<Vec4i> lines;
    Size size = src0.size();
    vector <double> angle_samples;

    // finds line segments in a binary image using the probabilistic Hough transform.
    // threshold – Accumulator threshold parameter. Only those lines are returned that get enough votes().
    // minLineLength – Minimum line length.Line segments shorter than that are rejected.
    // maxLineGap – Maximum allowed gap between points on the same line to link them.
    HoughLinesP(src0, lines, 1, CV_PI / 180, threshold, minLineLenght, maxLineGap);

    Mat disp_lines(size, CV_8UC1, Scalar(0, 0, 0));
    double angle = 0.;
    unsigned nb_lines = (int)lines.size();
    for (unsigned i = 0; i < nb_lines; ++i)
    {

        line(disp_lines, Point(lines[i][0], lines[i][1]),
             Point(lines[i][2], lines[i][3]), Scalar(255, 0, 0));
        double cur_angle = atan2((double)lines[i][3] - lines[i][1],
                                 (double)lines[i][2] - lines[i][0]);
        angle_samples.push_back(cur_angle);
        angle += cur_angle;

    }

    if (display_on)
    {
        imshow("houghlines", disp_lines);
        imwrite("houghlines.jpg", disp_lines);
        if (nb_lines == 0)
            cout << "Very Small label resolution, nb_lines = 0" << endl;
    }

    if (nb_lines == 0)
    {
        return(angle);
    }


    angle /= nb_lines; // mean angle, in radians.
    angle = makehist(angle_samples, 60, CV_PI, -CV_PI);
    angle = ((angle * 180 / CV_PI));

    if (angle > 0)
        angle = angle - 180;
    if (display_on)
    {
        cout << "ANGLE : " << angle << endl;
    }
    if ((angle > -min_angle_to_assign_zero) && (angle < min_angle_to_assign_zero))
        angle = 0;
    else if ((angle < -max_angle_to_assign_zero) || (angle > max_angle_to_assign_zero))
        angle = 0;

    return(angle);
}
