#include <blending/naive_blender.h>

using namespace blending;

NaiveBlender::NaiveBlender(const std::vector<cv::Mat>& masks) : masks(masks) {
}

void NaiveBlender::Blend(const std::vector<cv::Mat>& images, cv::Mat& result) {
    for (int iter = 0; iter < images.size(); iter++) {
        images[iter].copyTo(result, masks[iter]);
    }
}
