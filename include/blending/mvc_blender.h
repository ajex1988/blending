#ifndef INCLUDE_BLENDING_MVCBLENDER_H_
#define INCLUDE_BLENDING_MVCBLENDER_H_

#include <blending/blender.h>

#include <opencv2/opencv.hpp>
#include <vector>

namespace blending {

class MvcBlender : public Blender {
public:
    MvcBlender(const std::vector<cv::Mat>& masks, const std::vector<std::vector<cv::Point> >& boundary, const std::vector<std::vector<int> >& boundary_idx, int pano_w, int pano_h);
    ~MvcBlender();
    void Blend(const std::vector<cv::Mat>& images, cv::Mat& result);
    static double AngleBetweenVector(cv::Point v1, cv::Point v2);
    static cv::Vec3d TriangleInterpolation(cv::Point p, cv::Point p1, cv::Point p2, cv::Point p3);
protected:
    void Triangulate(const std::vector<cv::Point>& bound, std::vector<cv::Point>& pts, std::vector<int>& tris, const cv::Mat& mask);
    void ComputeCoords();
    void ComputeCoords(const std::vector<cv::Point>& vertex, const std::vector<cv::Point>& boundary, const std::vector<int>& seam_idx, std::vector<double>& coords);
    void ComputeTriangle();

private:
    int pano_w, pano_h;
    // pixel data
    cv::Mat triangle_map;
    cv::Mat triangle_component;
    // boundary data
    std::vector<std::vector<cv::Point> > boundaries; // boundary points
    std::vector<std::vector<int> > seam_elements; // MVC boundary index
    // triangle data
    std::vector<std::vector<cv::Point> > vertexes; // Delaunay triangle vertexes
    std::vector<std::vector<cv::Point> > diff_vertexes; // MVC diff points
    std::vector<std::vector<int> > triangle_elements; // Triangle vertex index
    // mvc data
    std::vector<std::vector<double> > mvc_coords;
    std::vector<std::vector<double> > mvc_diff_coords;
};

}


#endif //INCLUDE_BLENDING_MVCBLENDER_H_
