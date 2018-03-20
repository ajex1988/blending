#ifndef INCLUDE_BLENDING_FEATHERBLENDER_H_
#define INCLUDE_BLENDING_FEATHERBLENDER_H_

#include <blending/blender.h>

#include <opencv2/opencv.hpp>
#include <vector>

namespace blending {

class FeatherBlender : public Blender {
public:
    FeatherBlender(const std::vector<cv::Mat>& masks, float sharpness = 0.02f);
    void Blend(const std::vector<cv::Mat> & image, cv::Mat& result);
private:
    int stitch_num;
    cv::Rect dst_roi;
    cv::Mat dst_weight;
    std::vector<cv::Mat> masks, weights;
};

}


#endif //INCLUDE_BLENDING_FEATHERBLENDER_H_
