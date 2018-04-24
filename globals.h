#pragma once
// Rectangles associated with circles
#define RECT_DEF_TOL 25
#define RECT_REDUCTION 0.25

// Probabilistic Hough Lines
#define MIN_NUM_VOTES_LINES 100
#define MIN_LINE_LENGTH 100.
#define MAX_LINE_GAP 100.
// Aspect Ratio of labels
#define LABEL_RATIO_MIN 2.0
#define LABEL_RATIO_MAX 2.6
#define LABEL_ANGLE_90RAD 0.075
//  Label Rotation
#define MIN_ANGLE_TO_ASSIGN_ZERO 5.
#define MAX_ANGLE_TO_ASSIGN_ZERO 175.
// Label Cropping
#define ROOM_FOR_ROTATION_IN_PIXELS 6.0
#define ROOM_FOR_CODE_EXTRACTION_IN_PIXELS 1.75
// Label Dimensions
#define LABEL_WIDTH_PIXELS 750
#define LABEL_HEIGHT_PIXELS 325
#define LABEL_VARIATION 100
// Hough Circle Constants
#define MIN_RADIUS 30
#define MAX_RADIUS 60
#define CIRCLE_DIST 125
#define MIN_NUM_VOTES_CIRCLES 100
#define MIN_NUM_VOTES_CIRCLES_IN_LABELS 80
// Threshold to start contrast stretching
#define STRETCH_THRESHOLD_IN_BITS 4
// Clean-up of undesirable blobs
#define MIN_BLOB_AREA 400
#define MAX_BLOB_AREA 1600
// Scales searches
#define LABEL_SCALES .10
// Number of times to search labels
#define LABEL_SEARCH_TRIALS 2
// Definition of night and day
#define BRIGHTNESS_LEVEL 0.5
// Tesseract programming
#define TESSERACT_PROG "bazaar"

typedef struct TAG
{
    int width_pxls;
    int height_pxls;
    int variation_pxls;
    int min_radius;
    int max_radius;
    int circle_dist;
    int min_blob_area;
    int max_blob_area;
    int stretch;
    double scale_search;
    int search_trials;
} Tag;
