#include "geometry.h"
#include "blobcleanup.h"


int ConvertQuadinRects(const Mat src, const vector<Vec3f>& Quadrilateral, Rect & Rectangle)
{

    Mat image = src.clone();
    vector<Vec3f> circles;
    vector<Vec3f> circles_quad;
    vector<Point> Quad;
    int X, Y;
    int maxX = 0;
    int minX = numeric_limits<int>::max();
    int maxY = 0;
    int minY = numeric_limits<int>::max();
    int actual_radius = 0;
    for (size_t j = 0; j < Quadrilateral.size(); j++)
    {
        Vec3f point = Quadrilateral[j];
        X = (int)point[0];
        Y = (int)point[1];
        actual_radius += (int)point[2];
        if (X > maxX)
        {
            maxX = X;
        }
        if (X < minX)
        {
            minX = X;
        }
        if (Y > maxY)
        {
            maxY = Y;
        }
        if (Y < minY)
        {
            minY = Y;
        }
    }
    actual_radius /= 4;
    Rect rect;

    rect.x = minX - (int)(ROOM_FOR_ROTATION_IN_PIXELS *actual_radius);
    rect.y = minY - (int)(ROOM_FOR_ROTATION_IN_PIXELS *actual_radius);
    rect.width  = (maxX - minX) + (int)(2*ROOM_FOR_ROTATION_IN_PIXELS *actual_radius);
    rect.height = (maxY - minY) + (int)(2*ROOM_FOR_ROTATION_IN_PIXELS *actual_radius);

    if (rect.x < 0)
        rect.x = 0;
    if (rect.y < 0)
        rect.y = 0;

    if ((rect.x + rect.width) > src.cols)
    {
        rect.width = src.cols - rect.x;
    }

    if ((rect.y + rect.height) > src.rows)
    {
        rect.height = src.rows - rect.y;
    }

    Rectangle = rect;
    return(actual_radius);
}
int ConvertQuadinLabelBounds(const Mat src, const vector<Vec3f>& Quadrilateral, Rect & Rectangle)
{
    int X, Y;
    int maxX = 0;
    int minX = numeric_limits<int>::max();
    int maxY = 0;
    int minY = numeric_limits<int>::max();
    int actual_radius = 0;
    for (size_t j = 0; j < Quadrilateral.size(); j++)
    {
        Vec3f point = Quadrilateral[j];
        X = (int)point[0];
        Y = (int)point[1];
        actual_radius += (int)point[2];
        if (X > maxX)
        {
            maxX = X;
        }
        if (X < minX)
        {
            minX = X;
        }
        if (Y > maxY)
        {
            maxY = Y;
        }
        if (Y < minY)
        {
            minY = Y;
        }
    }
    actual_radius /= 4;
    Rect rect;

    rect.x = minX ;
    rect.y = minY + (int) (round(ROOM_FOR_CODE_EXTRACTION_IN_PIXELS*actual_radius));
    rect.width  = (maxX - minX) ;
    rect.height = (maxY - minY) - (int) (2. * round(ROOM_FOR_CODE_EXTRACTION_IN_PIXELS*actual_radius));

    if (rect.x < 0)
        rect.x = 0;
    if (rect.y < 0)
        rect.y = 0;
    if (rect.x > src.cols)
        rect.x = 0;
    if (rect.y > src.rows)
        rect.y = 0;

    if ((rect.x + rect.width) > src.cols)
    {
        rect.width = src.cols - rect.x;
    }

    if ((rect.y + rect.height) > src.rows)
    {
        rect.height = src.rows - rect.y;
    }

    Rectangle = rect;
    return(actual_radius);
}


Mat Rotate(Mat src, double angle)
{

    Mat dst;
    Point2f pt(src.cols / 2.f, src.rows / 2.f);
    Mat r = getRotationMatrix2D(pt, angle, 1.0);
    warpAffine(src, dst, r, Size(src.cols,src.rows), cv::INTER_LANCZOS4);
    return dst;

}

void insertionSort2(double mag[], double ang[], int size)
{
    int i, j;
    double tmp;

    for (int k = 0; k < size; k++)
    {
        i = 0;
        j = i + 1;
        while (i < size - 1)
        {
            while (j < size)
            {
                if (mag[i] > mag[j])
                {
                    tmp = mag[i];
                    mag[i] = mag[j];
                    mag[j] = tmp;
                    tmp = ang[i];
                    ang[i] = ang[j];
                    ang[j] = tmp;
                }
                j++;
            }
            i++;
            j = i + 1;
        }
    }

}
void insertionSort1(double mag[], int size)
{
    int i, j;
    double tmp;

    for (int k = 0; k < size; k++)
    {
        i = 0;
        j = i + 1;
        while (i < size - 1)
        {
            while (j < size)
            {
                if (mag[i] > mag[j])
                {
                    tmp = mag[i];
                    mag[i] = mag[j];
                    mag[j] = tmp;
                }
                j++;
            }
            i++;
            j = i + 1;
        }
    }

}
double CalcAngle(Mat src_gray, const vector<Vec3f>& Quadrilateral, double min_angle_to_assign_zero, double max_angle_to_assign_zero)
{
    vector<Vec3f> vertices;
    double angles[6],mag[6];
    Vec3f point1, point2;
    double X1, Y1, X2, Y2, dx, dy;


    for (int j = 0; j < (int)Quadrilateral.size(); j++)
    {
        Vec3f point = Quadrilateral[j];
        vertices.push_back(point);
    }

    vertices.push_back(Quadrilateral[0]);

    for (int j = 0; j < (int)vertices.size()-1; j++)
    {
        point1 = vertices[j];
        point2 = vertices[j+1];
        X1 = point1[0];
        Y1 = point1[1];
        X2 = point2[0];
        Y2 = point2[1];
        dx = X1 - X2;
        dy = Y1 - Y2;
        angles[j] = ((atan2(dy, dx) * 180 / CV_PI));
        mag[j] = sqrt(dx * dx + dy * dy);
    }
    point1 = vertices[3];
    point2 = vertices[1];
    X1 = point1[0];
    Y1 = point1[1];
    X2 = point2[0];
    Y2 = point2[1];
    dx = X1 - X2;
    dy = Y1 - Y2;

    angles[4] = ((atan2(dy, dx) * 180 / CV_PI));
    mag[4] = sqrt(dx * dx + dy * dy);

    point1 = vertices[0];
    point2 = vertices[2];
    X1 = point1[0];
    Y1 = point1[1];
    X2 = point2[0];
    Y2 = point2[1];
    dx = X1 - X2;
    dy = Y1 - Y2;

    angles[5] = ((atan2(dy, dx) * 180 / CV_PI));
    mag[5] = sqrt(dx * dx + dy * dy);
    insertionSort2(mag, angles, 6);

    // take the median value
    double angle = angles[3];

    if (angle > 0)
    {
        angle = angle - 180;
    }

    if ((angle > -min_angle_to_assign_zero) && (angle < min_angle_to_assign_zero))
        angle = 0;
    else if ((angle < -max_angle_to_assign_zero) || (angle > max_angle_to_assign_zero))
        angle = 0;

    return(angle);


}


Rect ReduceRect(Rect R, double val)
{
    Rect Rectangle;

    Rectangle.x = R.x + (int)(val * R.width);
    Rectangle.y = R.y + (int)(val * R.height);
    Rectangle.width = (int)(2.* val * R.width);
    Rectangle.height = (int)(2 * val * R.height);
    return(Rectangle);
}

Rect IncreaseRect(Mat src, Rect R, double val)
{
    int x, y, width, height;

    x = R.x + (int)(val * R.width);
    y = R.y + (int)(val * R.height);
    width = (int)((1. - val) * R.width);
    height = (int)((1. - val) * R.height);
    if (x < 0)
        x = 0;
    if (y < 0)
        y = 0;
    if ((x + width) > src.cols)
        width = src.cols - x;
    if ((y + height) > src.rows)
        height = src.rows - y;

    Rect Rectangle = { x,y,width,height };

    return(Rectangle);
}



int CircularDetection(Mat &img, bool contrast_stretch, const int min_dist_center, const int min_number_of_votes, const int min_radius, const int max_radius, vector<Vec3f> &circles)
{
    Mat dst = img.clone();


    if (contrast_stretch==false)
    {
        GaussianBlur(img, dst, cv::Size(9, 9), 2.5, 2.5);
    }
#ifdef __GNUC__
    else
    {
        GaussianBlur(img, dst, cv::Size(7, 7), 0.55, 0.55);
    }
#endif
    HoughCircles(dst, circles, CV_HOUGH_GRADIENT,
                 2,   // accumulator resolution (size of the image / 2)
                 min_dist_center,  // minimum distance between two circles
                 20, // Canny high threshold
                 min_number_of_votes, // minimum number of votes
                 min_radius, max_radius); // min and max radius

    return((int)circles.size());

}

int CircularDetection(Mat &img, const int min_dist_center, const int min_number_of_votes, const int min_radius, const int max_radius, vector<Vec3f> &circles)
{
    Mat dst = img.clone();


    GaussianBlur(img, dst, cv::Size(7, 7), 1.5, 1.5);


    HoughCircles(dst, circles, CV_HOUGH_GRADIENT,
                 2,   // accumulator resolution (size of the image / 2)
                 min_dist_center,  // minimum distance between two circles
                 20., // Canny high threshold
                 min_number_of_votes, // minimum number of votes
                 min_radius, max_radius); // min and max radius

    return((int)circles.size());

}

void DrawCircles(const Mat image, vector<Vec3f> &circles, string name, double scale)
{
    Mat src = image.clone();
    const char* window_name = name.c_str();

    std::vector<cv::Vec3f>::
    const_iterator itc = circles.begin();


    while (itc != circles.end())
    {

        circle(src,
               Point((int)(*itc)[0], (int)(*itc)[1]), // circle centre
               (int)(*itc)[2],                            // circle radius
               Scalar(0, 255, 0),                           // color
               5);                                        // thickness

        ++itc;
    }
    namedWindow(window_name, 1);
    cv::namedWindow(window_name, WINDOW_NORMAL);
    // resize to fit display
    Mat resized;
    int sx = (int)(scale*image.cols);
    int sy = (int)(scale*image.rows);
    resize(src, resized, Size(sx,sy));
    imshow(window_name, resized);
}

void ShowFrame(const Mat image, string name, double scale)
{
    Mat src = image.clone();
    const char* window_name = name.c_str();


    namedWindow(window_name, 1);
    cv::namedWindow(window_name, WINDOW_NORMAL);
    // resize to fit display
    Mat resized;
    int sx = (int)(scale*image.cols);
    int sy = (int)(scale*image.rows);
    resize(src, resized, Size(sx, sy));
    imshow(window_name, resized);
}

void PrintCircles(vector<Vec3f> &circles)
{

    std::vector<cv::Vec3f>::
    const_iterator itc = circles.begin();


    while (itc != circles.end())
    {
        cout << "X=" << (int)(*itc)[0] << ',';
        cout << "Y=" << (int)(*itc)[1] << ',' << endl;
        ++itc;
    }
}

// finds a cosine of angle between vectors
// from pt0->pt1 and from pt0->pt2
static double angle(Vec3f pt1, Vec3f pt2, Vec3f pt0)
{
    double dx1 = pt1[0] - pt0[0];
    double dy1 = pt1[1] - pt0[1];
    double dx2 = pt2[0] - pt0[0];
    double dy2 = pt2[1] - pt0[1];
    return ((dx1*dx2 + dy1*dy2) / sqrt((dx1*dx1 + dy1*dy1)*(dx2*dx2 + dy2*dy2) + 1e-10));
}
static double dist(Vec3f point1, Vec3f point2)
{
    double d = sqrt((point1[0] - point2[0]) *  (point1[0] - point2[0]) +
                    (point1[1] - point2[1]) *  (point1[1] - point2[1]));

    return(d);
}

void comb_with_repetions(int N, int K, vector<vector<int>> &indexset)
{
    std::string bitmask(K, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's

    do
    {
        vector<int> buffer;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
            {
                //std::cout << " " << i;
                buffer.push_back(i);
            }
        }
        // find possible arrangements
        vector<int> arrange = buffer;
        sort(arrange.begin(), arrange.end());
        do
        {
            indexset.push_back(arrange);
        }
        while (std::next_permutation(arrange.begin(), arrange.end()));

        //std::cout << std::endl;
    }
    while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}
bool FindQuad4pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad, int &scale)
{
    if (circles.size() < 4)
        return(false);
    int N = (int)circles.size();
    int K = 4;

    vector<Vec3f> Buffer(4);
    bool cond1 = false;
    bool cond2 = false;
    bool cond3 = false;
    bool cond4 = false;
    vector<vector<int>> indexset;
    int count, number_of_trials = tag.search_trials;
    int width_pxls = tag.width_pxls;
    int height_pxls = tag.height_pxls;
    int variation_pxls = tag.variation_pxls;
    double rects_scales = tag.scale_search;
    // Generate Combinations with repetitions table
    comb_with_repetions(N,K,indexset);
    scale = -1;

    count = 0;
    do
    {
        for (size_t i = 0; i < indexset.size(); i++)
        {
            vector<int> cset = indexset[i];
            Vec3f pointa = circles[cset[0]];
            Vec3f pointb = circles[cset[1]];
            Vec3f pointc = circles[cset[2]];
            Vec3f pointd = circles[cset[3]];
            if (fabs(dist(pointa, pointb) - dist(pointc, pointd)) < tol)
                cond1 = true;
            if (fabs(dist(pointa, pointc) - dist(pointb, pointd)) < tol)
                cond2 = true;
            if (fabs(dist(pointa, pointd) - dist(pointb, pointc)) < tol)
                cond3 = true;

            Buffer[0] = pointa;
            Buffer[1] = pointb;
            Buffer[2] = pointc;
            Buffer[3] = pointd;

            cond4 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);

            if ((cond1 == true) && (cond2 == true) && (cond3 == true) && (cond4 == true))
            {
                Quad.push_back(pointa);
                Quad.push_back(pointb);
                Quad.push_back(pointc);
                Quad.push_back(pointd);
                scale = count;
                return(true);
            }
            // reset conditions
            cond1 = false;
            cond2 = false;
            cond3 = false;
            cond4 = false;
        }
        // readjust search the tag may have smaller dimensions
        count++;
        width_pxls -= (int)(rects_scales*tag.width_pxls);
        height_pxls -= (int)(rects_scales*tag.height_pxls);
        variation_pxls -= (int)(rects_scales*tag.variation_pxls);

        cond1 = false;
        cond2 = false;
        cond3 = false;
        cond4 = false;

        if ((width_pxls < 0) || (height_pxls < 0) || (variation_pxls < 0))
        {
            return(false);
        }
    }
    while (count < number_of_trials);

    return(false);
}


bool FindQuad3pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad, int &scale)
{

    if (circles.size() < 3)
        return(false);
    int N = (int)circles.size();
    int K = 3;
    vector<Vec3f> Buffer(4);
    bool cond1 = false;
    bool cond2 = false;
    bool cond3 = false;
    vector<vector<int>> indexset;
    int count, number_of_trials = tag.search_trials;
    int width_pxls = tag.width_pxls;
    int height_pxls = tag.height_pxls;
    int variation_pxls = tag.variation_pxls;
    double rects_scales = tag.scale_search;
    // Generate Combinations with repetitions table
    comb_with_repetions(N, K, indexset);
    scale = -1;

    //int xsize = (int) indexset.size();
    // cout << "xsize " << xsize << endl;
    count = 0;
    do
    {
        for (size_t i = 0; i < indexset.size(); i++)
        {
            vector<int> cset = indexset[i];
            Vec3f pointa = circles[cset[0]];
            Vec3f pointb = circles[cset[1]];
            Vec3f pointc = circles[cset[2]];
            double distab = dist(pointa, pointb);
            double distbc = dist(pointb, pointc);
            double distac = dist(pointa, pointc);

            double diagonal = MAX(distab, MAX(distbc, distac));
            double height = MIN(distab, MIN(distbc, distac));
            double width = sqrt(diagonal*diagonal - height*height);


            double rratio = width / height;
            double alpha = fabs(angle(pointa, pointc, pointb));

#if 0
            if (alpha < 0.1)
            {
                cout << "distab = " << distab << endl;
                cout << "distbc = " << distbc << endl;
                cout << "distac = " << distac << endl;
                cout << width << " " << height << " " << diagonal << endl;
                cout << "angle =" << alpha << endl;
            }
#endif

            if (alpha < LABEL_ANGLE_90RAD)
            {
                cond2 = true;
                if ((rratio >= LABEL_RATIO_MIN) && (rratio < LABEL_RATIO_MAX))
                    cond1 = true;
            }

            if ((cond1 == true) && (cond2 == true))
            {
                Vec3f diffab, diffbc, pointd;

                diffab[0] = pointb[0] - pointa[0];
                diffab[1] = pointb[1] - pointa[1];
                diffbc[0] = pointc[0] - pointb[0];
                diffbc[1] = pointc[1] - pointb[1];

                pointd = pointc - diffab;
                pointd[2] = (pointa[2] + pointb[2] + pointc[2]) / 3;

                Buffer[0] = pointa;
                Buffer[1] = pointb;
                Buffer[2] = pointc;
                Buffer[3] = pointd;
#if 0
                cout << "--------" << endl;
                cout << "A =" << pointa << endl;
                cout << "B =" << pointb << endl;
                cout << "C =" << pointc << endl;
                cout << "D =" << pointd << endl;
                cout << "--------" << endl;

#endif
                cond3 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);
                if (cond3 == true)
                {
                    Quad.push_back(pointa);
                    Quad.push_back(pointb);
                    Quad.push_back(pointc);
                    Quad.push_back(pointd);
                    scale = count;
                    return(true);
                }
                // this routine above checks if rect is horizontally aligned
                Vec3f diffba = -diffab;
                Vec3f diffbd = diffba + diffbc;
                pointd[0] = diffbd[0] + pointb[0];
                pointd[1] = diffbd[1] + pointb[1];
                pointd[2] = (pointa[2] + pointb[2] + pointc[2]) / 3;

                Buffer[0] = pointa;
                Buffer[1] = pointb;
                Buffer[2] = pointc;
                Buffer[3] = pointd;

                cond3 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);

                if (cond3 == true)
                {
                    Quad.push_back(pointa);
                    Quad.push_back(pointb);
                    Quad.push_back(pointc);
                    Quad.push_back(pointd);
                    scale = count;
                    return(true);
                }
            }
            // reset conditions
            cond1 = false;
            cond2 = false;
            cond3 = false;
        }
        // readjust search the tag may have smaller dimensions
        count++;
        width_pxls -= (int)(rects_scales*tag.width_pxls);
        height_pxls -= (int)(rects_scales*tag.height_pxls);
        variation_pxls -= (int)(rects_scales*tag.variation_pxls);

        cond1 = false;
        cond2 = false;
        cond3 = false;
        if ((width_pxls < 0) || (height_pxls < 0) || (variation_pxls < 0))
        {
            return(false);
        }
    }
    while (count < number_of_trials);
    return(false);
}

bool isValidQuad(vector<Vec3f> &Quad,double qwidth, double qheight, int tol)
{
    vector<Vec3f> vertices;
    double mag[6];
    Vec3f point1, point2;
    double X1, Y1, X2, Y2, dx, dy;
    bool cond = false;
    bool cond1 = false;
    bool cond2 = false;
    bool cond3 = false;

    if (Quad.size() != 4)
        return(false);
    for (int j = 0; j < (int)Quad.size(); j++)
    {
        Vec3f point = Quad[j];
        vertices.push_back(point);
    }

    vertices.push_back(Quad[0]);

    for (int j = 0; j < (int)vertices.size() - 1; j++)
    {
        point1 = vertices[j];
        point2 = vertices[j + 1];
        X1 = point1[0];
        Y1 = point1[1];
        X2 = point2[0];
        Y2 = point2[1];
        dx = X1 - X2;
        dy = Y1 - Y2;
        mag[j] = sqrt(dx * dx + dy * dy);
    }
    point1 = vertices[3];
    point2 = vertices[1];
    X1 = point1[0];
    Y1 = point1[1];
    X2 = point2[0];
    Y2 = point2[1];
    dx = X1 - X2;
    dy = Y1 - Y2;

    mag[4] = sqrt(dx * dx + dy * dy);

    point1 = vertices[0];
    point2 = vertices[2];
    X1 = point1[0];
    Y1 = point1[1];
    X2 = point2[0];
    Y2 = point2[1];
    dx = X1 - X2;
    dy = Y1 - Y2;

    mag[5] = sqrt(dx * dx + dy * dy);
    insertionSort1(mag,6);
    double height   = (mag[0]+mag[1])/2;
    double width    = (mag[2]+mag[3])/2;
    //double diagonal = 0.5 * (mag[4] + mag[5]);
#if 0
    cout << "Height: " << height << endl;
    cout << "width: " << width << endl;
    //cout << "diagonal: " << diagonal << endl;
#endif


    if (((height - qheight)*(height - qheight)) < tol*tol)
        cond1 = true;
    else
        return(false);
    if (((width - qwidth)*(width - qwidth)) < tol*tol)
        cond2 = true;
    else
        return(false);

    double maxCosine = 0.;
    for (int j = 2; j < 5; j++)
    {
        // find the maximum cosine of the angle between joint edges
        double cosine = fabs(angle(vertices[j % 4], vertices[j - 2], vertices[j - 1]));
        maxCosine = MAX(maxCosine, cosine);
    }

    // if cosines of all angles are small
    // (all angles are ~90 degree) then write quandrangular
    // vertices to resultant sequence LABEL_ANGLE_90RAD
    // if (maxCosine <  LABEL_ANGLE_90RAD)
    //
    if (maxCosine < LABEL_ANGLE_90RAD)
        cond3 = true;
    else
        return false;

    if ((cond1) && (cond2) && (cond3))
        cond = true;
    else
        cond = false;

    return(cond);
}
#if 0
bool FindQuad4pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad)
{
    int size = (int)circles.size();
    int i, j, k, l;
    vector<Vec3f> Buffer(4);
    bool cond1 = false;
    bool cond2 = false;
    bool cond3 = false;
    bool cond4 = false;
    int count, number_of_trials = tag.search_trials;
    int width_pxls = tag.width_pxls;
    int height_pxls = tag.height_pxls;
    int variation_pxls = tag.variation_pxls;

    count = 0;
    do
    {
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                for (k = 0; k < size; k++)
                {
                    for (l = 0; l < size; l++)
                    {
                        if ((i != j) && (i != k) && (i != l) && (j != k) && (j != l) && (k != l))
                        {
                            Vec3f pointa = circles[i];
                            Vec3f pointb = circles[j];
                            Vec3f pointc = circles[k];
                            Vec3f pointd = circles[l];
                            if (fabs(dist(pointa, pointb) - dist(pointc, pointd)) < tol)
                                cond1 = true;
                            if (fabs(dist(pointa, pointc) - dist(pointb, pointd)) <  tol)
                                cond2 = true;
                            if (fabs(dist(pointa, pointd) - dist(pointb, pointc)) <  tol)
                                cond3 = true;

                            Buffer[0] = pointa;
                            Buffer[1] = pointb;
                            Buffer[2] = pointc;
                            Buffer[3] = pointd;

                            cond4 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);

                            if ((cond1 == true) && (cond2 == true) && (cond3 == true) && (cond4 == true))
                            {
                                Quad.push_back(pointa);
                                Quad.push_back(pointb);
                                Quad.push_back(pointc);
                                Quad.push_back(pointd);
                                return(true);
                            }
                            else
                            {
                                cond1 = false;
                                cond2 = false;
                                cond3 = false;
                                cond4 = false;
                            }
                        }
                    }
                }
            }
        }
        // readjust search parameters
        width_pxls -= variation_pxls;
        height_pxls -= variation_pxls;
        count++;
    }
    while (count < number_of_trials);

    return(false);
}

bool FindQuad3pts(vector<Vec3f> &circles, double tol, Tag tag, vector<Vec3f> &Quad)
{
    int size = (int)circles.size();
    int i, j, k;
    vector<Vec3f> Buffer(4);
    bool cond1 = false;
    bool cond2 = false;
    bool cond3 = false;
    int count, number_of_trials = tag.search_trials;
    int width_pxls = tag.width_pxls;
    int height_pxls = tag.height_pxls;
    int variation_pxls = tag.variation_pxls;

    count = 0;
    do
    {
        for (i = 0; i < size; i++)
        {
            for (j = 0; j < size; j++)
            {
                for (k = 0; k < size; k++)
                {
                    if ((i != j) && (i != k) && (j != k))
                    {
                        Vec3f pointa = circles[i];
                        Vec3f pointb = circles[j];
                        Vec3f pointc = circles[k];
                        double distab = dist(pointa, pointb);
                        double distbc = dist(pointb, pointc);
                        double distac = dist(pointa, pointc);
                        double diagonal = MAX(distab, MAX(distbc, distac));
                        double height = MIN(distab, MIN(distbc, distac));
                        double width = sqrt(diagonal*diagonal - height*height);

#if 0
                        cout << width << " " << height << " " << diagonal << endl;
#endif
                        double ratio = width / height;
                        double alpha = fabs(angle(pointa, pointc, pointb));

                        if (alpha < LABEL_ANGLE_90RAD)
                        {
                            cond2 = true;
                            if ((ratio >= LABEL_RATIO_MIN) && (ratio < LABEL_RATIO_MAX))
                                cond1 = true;
                        }
                        if ((cond1 == true) && (cond2 == true))
                        {

                            Vec3f diffab, diffbc, pointd;

                            diffab[0] = pointb[0] - pointa[0];
                            diffab[1] = pointb[1] - pointa[1];
                            diffbc[0] = pointc[0] - pointb[0];
                            diffbc[1] = pointc[1] - pointb[1];

                            pointd = pointc - diffab;
                            pointd[2] = (pointa[2] + pointb[2] + pointc[2]) / 3;

                            Buffer[0] = pointa;
                            Buffer[1] = pointb;
                            Buffer[2] = pointc;
                            Buffer[3] = pointd;
#if 0
                            cout << "--------" << endl;
                            cout << "A =" << pointa << endl;
                            cout << "B =" << pointb << endl;
                            cout << "C =" << pointc << endl;
                            cout << "D =" << pointd << endl;
                            cout << "--------" << endl;

#endif
                            cond3 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);
                            if (cond3 == true)
                            {
                                Quad.push_back(pointa);
                                Quad.push_back(pointb);
                                Quad.push_back(pointc);
                                Quad.push_back(pointd);
                                return(true);
                            }
                            // this routine above checks if rect is horizontally aligned
                            Vec3f diffba = -diffab;
                            Vec3f diffbd = diffba + diffbc;
                            pointd[0] = diffbd[0] + pointb[0];
                            pointd[1] = diffbd[1] + pointb[1];
                            pointd[2] = (pointa[2] + pointb[2] + pointc[2]) / 3;

                            Buffer[0] = pointa;
                            Buffer[1] = pointb;
                            Buffer[2] = pointc;
                            Buffer[3] = pointd;

                            cond3 = isValidQuad(Buffer, width_pxls, height_pxls, variation_pxls);

                            if (cond3 == true)
                            {
                                Quad.push_back(pointa);
                                Quad.push_back(pointb);
                                Quad.push_back(pointc);
                                Quad.push_back(pointd);
                                return(true);
                            }
                        }
                        else
                        {
                            cond1 = false;
                            cond2 = false;
                            cond3 = false;
                        }
                    }
                }
            }
        }
        count++;
        width_pxls -= variation_pxls;
        height_pxls -= variation_pxls;
        count++;
    }
    while (count < number_of_trials);
    return(false);
}
#endif
