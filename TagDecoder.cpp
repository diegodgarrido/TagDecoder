#include "Tagdecoder.h"
#include "opencv2/core/version.hpp"


int64 work_begin = 0;
int64 work_end = 0;


static void workBegin()
{
    work_begin = getTickCount();
}

static void workEnd()
{
    work_end = getTickCount() - work_begin;
}

static double getTime()
{
    return work_end / ((double)getTickFrequency())* 1000.;
}

float GetCode(tesseract::ResultIterator* ri, tesseract::PageIteratorLevel level, string &output, vector<float> &confidences_per_char,vector<char> &charstring)
{

    char code[32] = {};
    float code_confidence[32];

    int count = 0;
    if (ri != 0)
    {
        do
        {
            const char* symbol = ri->GetUTF8Text(level);
            //float conf = ri->Confidence(level);
            if (symbol != 0 && strlen(symbol) != 0)
            {
                tesseract::ChoiceIterator ci(*ri);
                float max = numeric_limits<float>::min();
                do
                {
                    const char* choice = ci.GetUTF8Text();
                    if (ci.Confidence() > max)
                    {
                        code[count] = *choice;
                        code_confidence[count] = ci.Confidence();
                        max = ci.Confidence();
                    }
                }
                while (ci.Next());
            }
            delete[] symbol;
            count++;
        }
        while ((ri->Next(level)));
    }


    string code_string(code);
    output = code_string;

    float min = numeric_limits<float>::max();
    for (int i = 0; i < count; i++)
    {
        charstring.push_back(code[i]);
        confidences_per_char.push_back(code_confidence[i]);
        if (code_confidence[i] < min)
        {
            min = code_confidence[i];
        }
    }

    return(min);
}

/** @function main */
int main(int argc, char** argv)
{
    char dirname[200]; // = { "C:\\Work\\TagDecoder" };
    string dirstring;
    vector<string> image_names;
    Rect boxbound;
    double min_angle_to_assign_zero = MIN_ANGLE_TO_ASSIGN_ZERO;
    double max_angle_to_assign_zero = MAX_ANGLE_TO_ASSIGN_ZERO;
    bool display_on = false;
    bool write_on = true;
    bool enable_flipped_text = true;
    int correct_labels = 0;
    string tessdir;
    string outtext, goldtext;
    string tagdimensions;
    Tag tag, rtag;


    // character recognition vocabulary
    //const char  *voc = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    //const char  *voc = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    //const char  *voc = "ABCDEFGHIJKLMNPQRSTUVWXYZ0123456789";
    const char  *voc = "ABCDEFGHIJKLMNPQRTUVWXYZ0123456789";
#ifdef __GNUC__
    char  tessdata[100] = "/home/diego/work/TagDecoder/tessdata";
#else
    char tessdata[100] = "\\Work\\TagDecoder\\tessdata";
#endif


    const char* keys =
        "{ h help                     |                   | print help message   }"
        "{ D dirname                  |                   | directory name (full path) where multiple images are (batch processing)}"
        "{ t tessdir                  |                   | directory name (full path) where tesseract is stored}"
        "{ b tagdimensions            |                   | tag characteristics   }"
        "{ i input_image              |  image.jpg        | specify input full path image or *.jpg  }"
        "{ o outtext                  |                   | decoding text file }"
        "{ g goldtext                 |                   | gold text file }"
        "{ a min_angle_to_assign_zero |                   | smallest value in degrees for not applying any rotation to label }"
        "{ A max_angle_to_assign_zero |                   | biggest value in degrees for not applying any rotation to label }"
        "{ f enable_flipped_text      |                   | enable flipped text analysis }"
        "{ w write_on (1 or 0)        |                   | w=1 write the final binarized label       }"
        "{ d display_on (1 or 0)      |                   | d=1 display image       }";


    CommandLineParser cmd(argc, argv, keys);

    if (cmd.has("help"))
    {
        cout << "Usage: tagdecoder [options]" << endl;
        cout << "Available options:" << endl;
        cmd.printMessage();
        return EXIT_SUCCESS;
    }

    if (cmd.has("tessdir"))
    {
        string parse_check;
        tessdir = cmd.get<string>("t");
#ifdef __GNUC__
        strcpy(tessdata, tessdir.c_str());
#else
        strcpy_s(tessdata, tessdir.c_str());
#endif
    }

    // Tesseract init and programming
    tesseract::TessBaseAPI *OCRTes = new tesseract::TessBaseAPI();
    cout << "Tesseract-ocr version: " << OCRTes->Version() << endl;
#ifdef __GNUC__
    char dictionary[sizeof(TESSERACT_PROG) - 1];
    memcpy( dictionary, TESSERACT_PROG, sizeof(dictionary));
#else
    char* dictionary = TESSERACT_PROG;
#endif
    char* config_array[1] = { dictionary };
    char** tessconfig = &config_array[0];


    if (OCRTes->Init(tessdata, "eng", tesseract::OEM_DEFAULT, tessconfig, 1, NULL, NULL, false))
    {
        cout << "Could not initiate the Tesseract OCR, most likey datapath tessdata is not correctly set." << endl;
        return EXIT_FAILURE;
    }

    // treat the image as a single word line
    tesseract::PageSegMode pagesegmode = static_cast<tesseract::PageSegMode>(8);
    OCRTes->SetPageSegMode(pagesegmode);

    OCRTes->SetVariable("tessedit_char_whitelist",voc);

    std::cout << "OpenCV version: "
              << CV_MAJOR_VERSION << "."
              << CV_MINOR_VERSION << "."
              << CV_SUBMINOR_VERSION
              << std::endl;

    ofstream outputFile;

    if (cmd.has("outtext"))
    {
        string parse_check;
        outtext = cmd.get<string>("o");
        outputFile.open(outtext);
    }

    ifstream goldFile;

    if (cmd.has("goldtext"))
    {
        string parse_check;
        goldtext = cmd.get<string>("g");
        goldFile.open(goldtext);
    }

    string input_imageName = cmd.get<string>("i");

    if (cmd.has("dirname"))
    {
        string parse_check;
        parse_check = cmd.get<string>("D");
        dirstring=parse_check;
#ifdef __GNUC__
        strcpy(dirname, parse_check.c_str());
#else
        strcpy_s(dirname, parse_check.c_str());
#endif
    }


    bool exists = input_imageName.find("*.jpg") != string::npos;

    if (exists) // Multiple jpeg files
    {
        getFrameSuffixes(dirname, ".jpg", image_names);
        for (size_t i = 0; i < image_names.size(); i++)
        {
            string full_path(dirname);
            input_imageName = image_names.at(i);
#ifdef __GNUC__
            image_names.at(i) = full_path + '/' + input_imageName;
#else
            image_names.at(i) = full_path + '\\' + input_imageName;
#endif
        }
    }
    else // Single full path file
    {
        image_names.push_back(input_imageName);
    }

    if (cmd.has("tagdimensions"))
    {
        string alabel = cmd.get<string>("b");
        vector<string> bbs;
        istringstream split(alabel);
        for (string temp; getline(split, temp, ':'); bbs.push_back(temp));
        tag.width_pxls = stoi(bbs[0]);
        tag.height_pxls = stoi(bbs[1]);
        tag.variation_pxls = stoi(bbs[2]);
        tag.min_radius = stoi(bbs[3]);
        tag.max_radius= stoi(bbs[4]);
        tag.circle_dist= stoi(bbs[5]);
        tag.min_blob_area= stoi(bbs[6]);
        tag.max_blob_area= stoi(bbs[7]);
        tag.stretch = stoi(bbs[8]);
        tag.scale_search = stof(bbs[9]);
        tag.search_trials = stoi(bbs[10]);
        memcpy(&rtag,&tag,sizeof(Tag));
    }
    else
    {
        tag.width_pxls = LABEL_WIDTH_PIXELS;
        tag.height_pxls = LABEL_HEIGHT_PIXELS;
        tag.variation_pxls = LABEL_VARIATION;
        tag.min_radius = MIN_RADIUS;
        tag.max_radius= MAX_RADIUS;
        tag.circle_dist= CIRCLE_DIST;
        tag.min_blob_area= MIN_BLOB_AREA;
        tag.max_blob_area= MAX_BLOB_AREA;
        tag.stretch = STRETCH_THRESHOLD_IN_BITS;
        tag.scale_search =LABEL_SCALES;
        tag.search_trials = LABEL_SEARCH_TRIALS;
        memcpy(&rtag,&tag,sizeof(Tag));
    }


    if (cmd.has("min_angle_to_assign_zero"))
    {
        string alabel = cmd.get<string>("a");
        min_angle_to_assign_zero = (stof(alabel));
    }
    if (cmd.has("max_angle_to_assign_zero"))
    {
        string alabel = cmd.get<string>("A");
        max_angle_to_assign_zero = (stof(alabel));
    }
    if (cmd.has("display_on"))
    {
        string alabel = cmd.get<string>("d");
        if (alabel == "1")
            display_on = true;
    }
    if (cmd.has("enable_flipped_text"))
    {
        string alabel = cmd.get<string>("f");
        if (alabel == "1")
            enable_flipped_text = true;
    }
    if (cmd.has("write_on"))
    {
        string alabel = cmd.get<string>("w");
        if (alabel == "1")
            write_on = true;
    }
    double algo_time = 0.;
    workBegin();
    for (size_t i = 0; i < image_names.size(); i++)
    {
        Mat img, src, src_gray, src_str;
        double brightness, contrast, saturation;
        int scale_used, scale_used_er;
        vector<Vec3f> circles, circles_er;
        vector<Vec3f> Quad, Quad_er;
        Rect Rectangle, Rectangle_er;
        string output, output_flipped;
        float confidence, confidence_flipped;
        vector<float> confidences_per_char, confidences_per_char_flipped;
        vector<char> charstring, charstring_flipped;
        int blob_area = 0;
        bool search_cond = false, search_cond_er = false, contrast_stretch = false;

        input_imageName = image_names.at(i);
        imread(input_imageName).copyTo(img);
        // Check if image can be loaded
        if (img.empty())
        {
            cout << "Couldn't load " << input_imageName << endl;
            cmd.printMessage();
            return EXIT_FAILURE;
        }
        cout << "Processing: " << input_imageName << endl;


        // Initialize arguments for border boundary conditions
        int top, bottom, left, right;
        int borderType = BORDER_CONSTANT;
        Scalar value = Scalar(255., 255., 255.);
        top = (int)(0.20*img.rows);
        bottom = (int)(0.20*img.rows);
        left = (int)(0.20*img.cols);
        right = (int)(0.20*img.cols);

        copyMakeBorder(img, src, abs(top), abs(bottom), abs(left), abs(right), borderType, value);


        cvtColor(src, src_gray, CV_BGR2GRAY);
        getImgPerception(src_gray, brightness, contrast, saturation);
        src_str = src_gray;


        // Find Circles in the image non-stretch

        circles.clear();
        CircularDetection(src_str, contrast_stretch, rtag.circle_dist, MIN_NUM_VOTES_CIRCLES, rtag.min_radius, rtag.max_radius, circles);

        if (display_on)
        {
            PrintCircles(circles);
            DrawCircles(src, circles, "CirclesI", 0.125);
        }

        // Find proper rectangle
        Quad.clear();
        search_cond = FindQuad4pts(circles, RECT_DEF_TOL, rtag, Quad, scale_used);


        if (display_on)
            DrawCircles(src, Quad, "CirclesR", 0.125);

        if (search_cond == false)
        {
            Quad.clear();
            search_cond = FindQuad3pts(circles, RECT_DEF_TOL, rtag, Quad, scale_used);
            cout << "(1)Doing a 3 pt search! " << endl;
        }

        // Find Circles in the image stretch
        if (search_cond == false)
        {
            int nbits = rtag.stretch;
            int threshold = 255 - (int)(pow(2., nbits));
            src_str = ConstrastStretching(src_gray, threshold, 0, 255, 255);
            contrast_stretch = true;
            cout << "Using contrast stretch! " << endl;
            circles.clear();
            CircularDetection(src_str, contrast_stretch, rtag.circle_dist, MIN_NUM_VOTES_CIRCLES, rtag.min_radius, rtag.max_radius, circles);

            if (display_on)
            {
                PrintCircles(circles);
                DrawCircles(src, circles, "CirclesI", 0.125);
            }

            // Find proper rectangle
            Quad.clear();
            search_cond = FindQuad4pts(circles, RECT_DEF_TOL, rtag, Quad, scale_used);


            if (display_on)
                DrawCircles(src, Quad, "CirclesR", 0.125);

            if (search_cond == false)
            {
                Quad.clear();
                search_cond = FindQuad3pts(circles, RECT_DEF_TOL, rtag, Quad, scale_used);
                cout << "(2)Doing a 3 pt search! " << endl;
            }


        }



        if (search_cond == false)
        {
            cout << "Scale (1st) :" << scale_used << endl;
            TagMessage("Could not locate label!", dirstring, input_imageName, (int)i, outputFile, outtext, goldFile, goldtext, rtag, output, confidence);
            continue;
        }
        else
        {
            // Try
            if (display_on)
            {
                cout << "xxxxxxxxxxxxx" << endl;
                cout << "Rectangle associated with label:" << endl;
                PrintCircles(Quad);
                DrawCircles(src, Quad, "CirclesF", 0.125);
            }
        }

        ConvertQuadinRects(src_gray, Quad, Rectangle);
        double angle = CalcAngle(src_gray, Quad, min_angle_to_assign_zero, max_angle_to_assign_zero);

        if (display_on)
        {
            cout << "Text skew (degrees) : " << angle << endl;
        }

        // Extract Label
        Mat label_roi(src_str, Rectangle);
        Mat label_rot = Rotate(label_roi, angle);

        if (display_on)
        {
            ShowFrame(label_roi, "Label Roi", .25);
            ShowFrame(label_rot, "Label Rotated", .25);

            imwrite("label_roi.jpg", label_roi);
            imwrite("label_rot.jpg", label_rot);
        }

        Mat label_er_gray = label_rot.clone();

        circles_er.clear();
        CircularDetection(label_er_gray, contrast_stretch, rtag.circle_dist, MIN_NUM_VOTES_CIRCLES_IN_LABELS, rtag.min_radius, rtag.max_radius, circles_er);


        if (display_on)
        {
            PrintCircles(circles_er);
            DrawCircles(label_er_gray, circles_er, "CirclesX", 0.25);
        }


        search_cond_er = FindQuad4pts(circles_er, RECT_DEF_TOL, rtag, Quad_er, scale_used_er);

        if (search_cond_er == false)
        {
            Quad_er.clear();
            search_cond_er = FindQuad3pts(circles_er, RECT_DEF_TOL, rtag, Quad_er, scale_used_er);
            if (search_cond_er == false)
            {
                TagMessage("Could not locate ER!", dirstring, input_imageName, (int)i, outputFile, outtext, goldFile, goldtext, rtag, output, confidence);
                continue;
            }

        }
        ConvertQuadinLabelBounds(label_rot, Quad_er, Rectangle_er);

        if (display_on)
            DrawCircles(label_rot, Quad_er, "CLabel", 0.25);

        // Produce the final binary label
        int iangle = (int)angle;

        blob_area = ((iangle == 0) ? rtag.min_blob_area : rtag.max_blob_area);

        Mat label(label_rot, Rectangle_er);


        Mat labelg, labele;

        // Apply Gamma
        if ((brightness > BRIGHTNESS_LEVEL) && (brightness <= 1.0))
        {
            if (contrast_stretch == false)
            {
                GaussianBlur(label, labelg, cv::Size(11, 11), 2.25, 2.25);
                labele = GammaCorrection(labelg, 0.45);
            }
            else
            {
                GaussianBlur(label, labelg, cv::Size(11, 11), 5.25, 5.25);
                labele = labelg;
            }
        }
        else
        {
            GaussianBlur(label, labelg, cv::Size(29, 29), 10, 10);
            labele = labelg; // GammaCorrection(labelg, 1.4);
        }
        label = labele;
        // Start binarization
        Mat mask = binarizebw(label, brightness, false);
        Mat final_mask = blobcleanup(mask, blob_area);


        if (display_on)
        {
            ShowFrame(label, "Label", .25);
            imshow("Binarized Label", mask);
            imshow("Binarized Label Clean", final_mask);
            imwrite("label.jpg", label);
        }
        // Run tesseract

        Rect textROI(0, 0, final_mask.cols, final_mask.rows);
        OCRTes->TesseractRect(final_mask.data, 1, (int)final_mask.step1(), textROI.x, textROI.y, textROI.width, textROI.height);
        confidence = GetCode(OCRTes->GetIterator(), tesseract::RIL_SYMBOL, output, confidences_per_char, charstring);

        if (display_on)
        {
            for (size_t cc = 0; cc < charstring.size(); cc++)
            {
                cout << "CHAR" << "[" << cc << "] =" << charstring[cc] << " -->Prob = " << confidences_per_char[cc] << endl;
            }
        }

        if (enable_flipped_text)
        {
            // Flip (Mirror) vertically;
            Mat final_mask_flipped = Mat(final_mask.rows, final_mask.cols, CV_8UC1, cvScalar(0));

            flip(final_mask, final_mask_flipped, -1);
            if (display_on)
            {
                imshow("Flipped Binarized Label", final_mask_flipped);
            }

            // Run tesseract on the flipped label
            Rect textROI(0, 0, final_mask_flipped.cols, final_mask_flipped.rows);
            OCRTes->TesseractRect(final_mask_flipped.data, 1, (int)final_mask_flipped.step1(), textROI.x, textROI.y, textROI.width, textROI.height);
            confidence_flipped = GetCode(OCRTes->GetIterator(), tesseract::RIL_SYMBOL, output_flipped, confidences_per_char_flipped, charstring_flipped);

            // Check confidences
            if ((output != "") && (output_flipped != ""))
            {
                if (confidence_flipped > confidence)
                {
                    final_mask = final_mask_flipped;
                    output = output_flipped;
                    confidence = confidence_flipped;
                }
            }
            else if ((output == "") && (output_flipped != ""))
            {
                final_mask = final_mask_flipped;
                output = output_flipped;
                confidence = confidence_flipped;
            }
            else
            {
                output = "empty";
                confidence = 0;
            }
            if (display_on)
            {
                for (size_t cc = 0; cc < charstring_flipped.size(); cc++)
                {
                    cout << "CHAR" << "[" << cc << "] =" << charstring_flipped[cc] << " -->Prob = " << confidences_per_char_flipped[cc] << endl;
                }
            }
        }
        cout << "OCR_Tesseract  output = " << output << endl;
        cout << "Confidence :" << confidence << endl;
        cout << "Scale (1st) :" << scale_used << endl;
        cout << "Scale (2nd) :" << scale_used_er << endl;

        if (display_on)
        {
            imshow("Final Binarized", final_mask);
        }

        if (write_on)
        {
            WriteMessage(final_mask, correct_labels, dirstring, input_imageName, (int)i, outputFile, outtext, goldFile, goldtext, rtag, output, confidence);
        }


    } // end of the loop
    if (outtext != "")
    {
        if (!goldtext.empty())
        {
            outputFile <<  "Correct Labels:  " << correct_labels << " " << 100. * ((float)correct_labels / (float)image_names.size()) << endl;
        }
        outputFile.close();
    }
    if (goldtext != "")
    {
        goldFile.close();
    }
    workEnd();
    algo_time = getTime();
    OCRTes->Clear();
    OCRTes->End();
    cout << "Label Finder run time: " << algo_time << " ms" << endl << "\n";
    cout << "Average Label Finder run time: " << algo_time / ((int)image_names.size() + 1) << " ms" << endl << "\n";

    char c = (char)waitKey(0);
    cout << c << endl;
    return EXIT_SUCCESS;
}
