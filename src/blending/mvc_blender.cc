#include <blending/mvc_blender.h>

#include <opencv2/opencv.hpp>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

#define EPS 1e-6

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CGAL::Delaunay_mesher_2<CDT, Criteria> Mesher;

typedef CDT::Vertex_handle Vertex_handle;
typedef K::Point_2   CGPoint;

using namespace blending;

void MvcBlender::Triangulate(const std::vector<cv::Point>& bound, std::vector<cv::Point>& pts, std::vector<int>& tris, const cv::Mat& mask) {
    std::vector<CGPoint> vertices;
    for (int i = 0; i < bound.size(); i++)
        vertices.push_back(CGPoint(bound[i].x, bound[i].y));
    CDT cdt;

    Vertex_handle v_start = cdt.insert(vertices[0]);
    Vertex_handle v_last = v_start;
    for (int i = 0; i < vertices.size(); i++) {
        Vertex_handle v_cur;
        if (i == vertices.size() - 1)
            v_cur = v_start;
        else
            v_cur = cdt.insert(vertices[i + 1]);

        cdt.insert_constraint(v_last, v_cur);
        v_last = v_cur;
    }

    assert(cdt.is_valid());

    Mesher mesher(cdt);
    mesher.refine_mesh();

    CDT::Finite_vertices_iterator vit;
    pts.clear(); tris.clear();
    for (vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); vit++) {
        pts.push_back(cv::Point(vit->point().hx(), vit->point().hy()));
    }

    CDT::Finite_faces_iterator fit;
    for (fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); fit++) {
        double cx = 0.0, cy = 0.0;
        bool out = false;
        for (int k = 0; k < 3; k++) {
            cv::Point pt(fit->vertex(k)->point().hx(), fit->vertex(k)->point().hy());
            cx += pt.x, cy += pt.y;
            if (mask.at<uchar>(pt.y, pt.x) == 0) {
                out = true;
                break;
            }
        }
        cx /= 3.0, cy /= 3.0;
        if (out || mask.at<uchar>(cvRound(cy), cvRound(cx)) == 0) continue;
        for (int k = 0; k < 3; k++) {
            cv::Point pt(fit->vertex(k)->point().hx(), fit->vertex(k)->point().hy());
            int index = find(pts.begin(), pts.end(), pt) - pts.begin();
            tris.push_back(index);
        }
    }
}

MvcBlender::MvcBlender(const std::vector<cv::Mat> & masks, const std::vector<std::vector<cv::Point> >& boundaries, const std::vector<std::vector<int> >& seams_idx, int pano_w, int pano_h)
    : boundaries(boundaries), seam_elements(seams_idx), pano_w(pano_w), pano_h(pano_h) {
    // allocate pixel data
    int blend_num = boundaries.size();
    vertexes.resize(blend_num);
    triangle_elements.resize(blend_num);
    mvc_coords.resize(blend_num);
    diff_vertexes.resize(blend_num - 1);
    mvc_diff_coords.resize(blend_num - 1);
    for (int iter = 0; iter < boundaries.size(); iter++) {
        const std::vector<cv::Point>& boundary = boundaries[iter];
        const std::vector<int>& seam_idx = seams_idx[iter];
        std::vector<cv::Point>& vertex = vertexes[iter];
        std::vector<int>& triangle = triangle_elements[iter];
        std::vector<double>& coords = mvc_coords[iter];

        Triangulate(boundary, vertex, triangle, masks[iter]);
        ComputeCoords(vertex, boundary, seam_idx, coords);
        if (iter < blend_num - 1) {
            std::vector<double>& diff_coords = mvc_diff_coords[iter];
            std::vector<cv::Point>& diff_vertex = diff_vertexes[iter];
            for (int k = iter + 1; k < blend_num; k++) {
                for (int i = 0; i < seams_idx[k].size(); i++) {
                    cv::Point p = boundaries[k][seams_idx[k][i]];
                    if (masks[iter].at<uchar>(p) == 0 && masks[iter].at<uchar>(p.y, (p.x + pano_w) % (2 * pano_w)) == 0) {
                        continue;
                    }
                    diff_vertex.push_back(p);
                }
            }

            ComputeCoords(diff_vertex, boundary, seam_idx, diff_coords);
        }
    }
    ComputeTriangle();
}

MvcBlender::~MvcBlender() {
}

void MvcBlender::Blend(const std::vector<cv::Mat>& images, cv::Mat& result_image) {
    int blend_num = boundaries.size();
    // precompute
    result_image = images[0].clone();
    std::vector<std::vector<cv::Vec3d> > membranes(blend_num);
    for (int iter = 0; iter < blend_num; iter++) {
        std::vector<cv::Vec3d> boundary_diff;
        for (int i = 0; i < seam_elements[iter].size(); i++) {
            cv::Point p = boundaries[iter][seam_elements[iter][i]];
            cv::Vec3b color_src = images[iter + 1].at<cv::Vec3b>(p.y, p.x % pano_w);
            cv::Vec3b color_dst = result_image.at<cv::Vec3b>(p.y, p.x % pano_w);
            cv::Vec3d diff = cv::Vec3d(color_dst[0] - color_src[0],
                               color_dst[1] - color_src[1],
                               color_dst[2] - color_src[2]);
            boundary_diff.push_back(diff);
        }
        std::vector<cv::Vec3d>& membrane = membranes[iter];
        int seam_size = seam_elements[iter].size();

        for (int i = 0; i < vertexes[iter].size(); i++) {
            cv::Vec3d offset(0.0, 0.0, 0.0);
            for (int j = 0; j < seam_size; j++) {
                offset += mvc_coords[iter][i * seam_size + j] * boundary_diff[j];
            }
            membrane.push_back(offset);
        }
        if (iter < blend_num - 1) {
            for (int i = 0; i < diff_vertexes[iter].size(); i++) {
                cv::Point p = diff_vertexes[iter][i];
                cv::Vec3d offset(0.0, 0.0, 0.0);
                for (int j = 0; j < seam_size; j++) {
                    offset += mvc_diff_coords[iter][i * seam_size + j] * boundary_diff[j];
                }
                cv::Vec3b color_src = images[iter + 1].at<cv::Vec3b>(p.y, p.x % pano_w);
                cv::Vec3b result_color;
                for (int k = 0; k < 3; k++) {
                    result_color[k] = uchar(MAX(MIN((color_src[k] + offset[k]), 255.0), 0.0));
                }
                result_image.at<cv::Vec3b>(p.y, p.x % pano_w) = result_color;
            }
        }
    }
    // parallel computing
    for (int y = 0; y < pano_h; y++) {
        for (int x = 0; x < pano_w; x++) {
            cv::Vec2i info = triangle_map.at<cv::Vec2i>(y, x);
            cv::Vec3d v = triangle_component.at<cv::Vec3d>(y, x);
            int image_idx = info[0];
            if (image_idx == 0) { // first image
                continue;
            }
            int triangle_idx = info[1];
            cv::Vec3d p1 = membranes[image_idx - 1][triangle_elements[image_idx - 1][triangle_idx]];
            cv::Vec3d p2 = membranes[image_idx - 1][triangle_elements[image_idx - 1][triangle_idx + 1]];
            cv::Vec3d p3 = membranes[image_idx - 1][triangle_elements[image_idx - 1][triangle_idx + 2]];
            cv::Vec3d color = p1 * v[0] + p2 * v[1] + p3 * v[2];

            cv::Vec3b color_src = images[image_idx].at<cv::Vec3b>(y, x);
            cv::Vec3b result_color;
            for (int k = 0; k < 3; k++) {
                result_color[k] = uchar(MAX(MIN((color_src[k] + color[k]), 255.0), 0.0));
            }
            result_image.at<cv::Vec3b>(y, x) = result_color;
            //result_image.at<cv::Vec3b>(y, x) = color_src;
        }
    }
}

double MvcBlender::AngleBetweenVector(cv::Point v1, cv::Point v2) {
    double l1 = sqrt(v1.ddot(v1)) + EPS, l2 = sqrt(v2.ddot(v2)) + EPS;
    double dot = v1.ddot(v2) / l1 / l2;
    return acos(MAX(MIN(dot, 1.0), -1.0));
}

cv::Vec3d MvcBlender::TriangleInterpolation(cv::Point p, cv::Point p1, cv::Point p2, cv::Point p3) {
    double s = (p2 - p1).cross(p3 - p1);
    double s1 = (p2 - p).cross(p3 - p) / s;
    double s2 = (p3 - p).cross(p1 - p) / s;
    double s3 = (p1 - p).cross(p2 - p) / s;
    return cv::Vec3d(s1, s2, s3);
}

void MvcBlender::ComputeCoords(const std::vector<cv::Point>& vertex, const std::vector<cv::Point>& boundary, const std::vector<int>& seam_idx, std::vector<double>& coords) {
    coords.resize(vertex.size() * seam_idx.size());
    int size = boundary.size();
    for (int i = 0; i < vertex.size(); i++) {
        double sum = 0.0;
        cv::Point px = vertex[i];
        for (int j = 0; j < seam_idx.size(); j ++) {
            int idx = seam_idx[j];
            cv::Point p0 = boundary[(idx + size - 1) % size];
            cv::Point p1 = boundary[idx];
            cv::Point p2 = boundary[(idx + size + 1) % size];
            double a1 = AngleBetweenVector(p0 - px, p1 - px);
            double a2 = AngleBetweenVector(p1 - px, p2 - px);
            double weight = (tan(abs(a1) / 2.0) + tan(abs(a2) / 2.0)) / (sqrt((p1 - px).dot(p1 - px)) + EPS);
            coords[i * seam_idx.size() + j] = weight;
            sum += weight;
        }
        assert(sum > 0.0);
        for (int j = 0; j < seam_idx.size(); j++) {
            coords[i * seam_idx.size() + j] /= sum;
        }
    }
}

void MvcBlender::ComputeTriangle() {
    triangle_map = cv::Mat::zeros(pano_h, pano_w, CV_32SC2);
    triangle_component = cv::Mat::zeros(pano_h, pano_w, CV_64FC3);
    for (int iter = 0; iter < boundaries.size(); iter++) {
        const std::vector<int> & triangle_list = triangle_elements[iter];
        const std::vector<cv::Point> & vertex = vertexes[iter];
        for (int i = 0; i < triangle_list.size(); i += 3) {
            cv::Point p1 = vertex[triangle_list[i]];
            cv::Point p2 = vertex[triangle_list[i + 1]];
            cv::Point p3 = vertex[triangle_list[i + 2]];

            int left = MIN(MIN(p1.x, p2.x), p3.x);
            int right = MAX(MAX(p1.x, p2.x), p3.x);
            int up = MIN(MIN(p1.y, p2.y), p3.y);
            int bottom = MAX(MAX(p1.y, p2.y), p3.y);

            std::vector<cv::Point> contour;
            contour.push_back(p1);
            contour.push_back(p2);
            contour.push_back(p3);
            contour.push_back(p1);

            for (int n = up; n <= bottom; ++n) {
                for (int m = left; m <= right; ++m) {
                    cv::Point p(m, n);
                    if (pointPolygonTest(contour, p, false) >= 0) {
                        triangle_map.at<cv::Vec2i>(n, m % pano_w) = cv::Vec2i(iter + 1, i);
                        triangle_component.at<cv::Vec3d>(n, m % pano_w) = TriangleInterpolation(p, p1, p2, p3);
                    }
                }
            }
        }
    }
}
