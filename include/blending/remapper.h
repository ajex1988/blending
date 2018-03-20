#ifndef INCLUDE_BLENDING_REMAPPER_H_
#define INCLUDE_BLENDING_REMAPPER_H_

#include <blending/script.h>

#include <vector>
#include <opencv2/opencv.hpp>

namespace blending {
class Remapper {
public:
    Remapper(std::string);
    ~Remapper();

    std::vector<cv::Mat> get_masks();

    void Remap(const std::vector<cv::Mat>& frames, std::vector<cv::Mat>& remap_frames);
private:
    /**
    * @brief initial the remap table for frame remapping
    */
    void initial_remap_table();
    void initial_transforms();

    void Transform(int idx, double, double, double*, double*);

    static double Cubic01(double x);
    static double Cubic12(double x);
private:
    int stitch_num;
    cv::Rect panorama_roi;
    Script ps;
    std::vector<TransformParam> transforms;
    std::vector<cv::Mat> masks;
    double ** remap_table_x, ** remap_table_y;
};

}
#endif // INCLUDE_BLENDING_REMAPPER_H_
