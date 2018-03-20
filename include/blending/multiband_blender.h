#ifndef INCLUDE_BLENDING_MULTIBANDBLENDER_H_
#define INCLUDE_BLENDING_MULTIBANDBLENDER_H_

#include <blending/blender.h>

#include <opencv2/opencv.hpp>
#include <vector>

namespace blending {

class MultiBandBlender : public Blender {
public:
    MultiBandBlender(const std::vector<cv::Mat>& masks);
    void Blend(const std::vector<cv::Mat> & images, cv::Mat& result);
protected:
    void BlendPyramdiLayers();
private:
    int stitch_num, band_num;
    cv::Rect dst_roi;
    cv::Mat dst_mask;
    std::vector<std::vector<cv::Mat> > src_laplace_pyr, src_weight_pyr;
    std::vector<cv::Mat> dst_laplace_pyr, dst_weight_pyr;
    std::vector<cv::Mat> masks;
};

}


#endif //INCLUDE_BLENDING_MULTIBANDBLENDER_H_
