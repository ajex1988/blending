#ifndef INCLUDE_BLENDING_BLENDER_H_
#define INCLUDE_BLENDING_BLENDER_H_

#include <opencv2/opencv.hpp>
#include <vector>

namespace blending {

class Blender {
public:
    virtual void Blend(const std::vector<cv::Mat> & image, cv::Mat& result) = 0;
};

}
#endif // INCLUDE_BLENDING_BLENDER_H_
