#ifndef INCLUDE_BLENDING_NAIVEBLENDER_H_
#define INCLUDE_BLENDING_NAIVEBLENDER_H_

#include <blending/blender.h>

#include <opencv2/opencv.hpp>
#include <vector>

namespace blending {

class NaiveBlender : public Blender {
public:
    NaiveBlender(const std::vector<cv::Mat>& masks);
    void Blend(const std::vector<cv::Mat> & image, cv::Mat& result);
private:
    std::vector<cv::Mat> masks;
};

}


#endif //INCLUDE_BLENDING_NAIVEBLENDER_H_
