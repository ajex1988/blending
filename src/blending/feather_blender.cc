#include <blending/feather_blender.h>

using namespace blending;

#define BLEND_WEIGHT_EPS 1e-6

FeatherBlender::FeatherBlender(const std::vector<cv::Mat>& masks, float sharpness) : masks(masks) {
    stitch_num = masks.size();
    for (auto& mask : masks) {
        dst_roi.width = std::max(mask.cols, dst_roi.width);
        dst_roi.height = std::max(mask.rows, dst_roi.height);
    }
    int w = dst_roi.width, h = dst_roi.height;
    weights.resize(stitch_num);
    dst_weight = cv::Mat::zeros(dst_roi.size(), CV_32F);
    for (int i = 0; i < stitch_num; i++) {
        distanceTransform(masks[i], weights[i], CV_DIST_L1, 3);
        threshold(weights[i] * sharpness, weights[i], 1.f, 1.f, cv::THRESH_TRUNC);
        for (int y = 0; y < h; ++y) {
            const float* weight_row = weights[i].ptr<float>(y);
            float* dst_weight_row = dst_weight.ptr<float>(y);
            for (int x = 0; x < w; ++x) {
                dst_weight_row[x] += weight_row[x];
            }
        }
    }
}

void FeatherBlender::Blend(const std::vector<cv::Mat>& images, cv::Mat& result) {
    int w = dst_roi.width, h = dst_roi.height;
    result = cv::Mat::zeros(h, w, CV_32FC3);
    for (size_t iter = 0; iter < stitch_num; iter ++) {
        for (int y = 0; y < h; ++y) {
            const cv::Point3_<uchar>* src_row = images[iter].ptr<cv::Point3_<uchar> >(y);
            const float* weight_row = weights[iter].ptr<float>(y);
            cv::Point3_<float>* dst_row = result.ptr<cv::Point3_<float> >(y);

            for (int x = 0; x < w; ++x) {
                dst_row[x].x += (float)src_row[x].x * weight_row[x];
                dst_row[x].y += (float)src_row[x].y * weight_row[x];
                dst_row[x].z += (float)src_row[x].z * weight_row[x];
            }
        }
    }
    for (int y = 0; y < h; ++y) {
        const float* dst_weight_row = dst_weight.ptr<float>(y);
        cv::Point3_<float>* dst_row = result.ptr<cv::Point3_<float> >(y);

        for (int x = 0; x < w; ++x) {
            dst_row[x].x /= dst_weight_row[x] + BLEND_WEIGHT_EPS;
            dst_row[x].y /= dst_weight_row[x] + BLEND_WEIGHT_EPS;
            dst_row[x].z /= dst_weight_row[x] + BLEND_WEIGHT_EPS;
        }
    }
    result.convertTo(result, CV_8U);
}
