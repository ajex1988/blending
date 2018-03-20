#include <blending/blending.h>
#include <utils/image_utils.h>

#define APP_NAME "mvc_blender_demo"

using namespace blending;

void CalculateBoundary(const std::vector<cv::Mat>& pano_masks, std::vector<cv::Mat>& blend_masks, std::vector<std::vector<cv::Point> >& boundaries, std::vector<std::vector<int> >& seams_idx) {
    int pano_h, pano_w;
    int stitch_num = pano_masks.size();
    if (!pano_masks.empty()) {
        pano_h = pano_masks[0].rows;
        pano_w = pano_masks[0].cols;
    }

    // copy to blend masks (size h * 2w)
    blend_masks.resize(stitch_num);
    for (int iter = 0; iter < stitch_num; iter++) {
        cv::Mat m = pano_masks[iter].clone();
        int w = m.cols, h = m.rows;
        blend_masks[iter] = cv::Mat::zeros(h, 2 * w, m.type());

        int ca = 0, cb = 0, ctop = 0;
        int xmin = m.cols, xmax = -1;
        for (int i = 0; i < m.rows; i++) {
            if (m.at<uchar>(i, 0) > 0) ca++;
            if (m.at<uchar>(i, m.cols - 1) > 0) cb++;
        }
        for (int i = 0; i < m.cols; i++) {
            if (m.at<uchar>(0, i) > 0) ctop++;
        }
        if (ca > 0 && cb > 0) {
            if (ctop == m.cols) { // sky
                pano_masks[iter].copyTo(blend_masks[iter](cv::Rect(0, 0, w, h)));
            }
            else { // cut
                cv::Rect left(0, 0, w / 2, h), right(w / 2, 0, w / 2, h);
                pano_masks[iter](right).copyTo(blend_masks[iter](right));
                pano_masks[iter](left).copyTo(blend_masks[iter](left + cv::Point(w, 0)));
            }
        }
        else {
            pano_masks[iter].copyTo(blend_masks[iter](cv::Rect(0, 0, w, h)));
        }
    }
    // dilate masks
    int dilation_size = 3;
    cv::Mat element = getStructuringElement(cv::MORPH_ELLIPSE,
        cv::Size(2 * dilation_size + 1, 2 * dilation_size + 1),
        cv::Point(dilation_size, dilation_size));
    for (int i = 0; i < stitch_num; i++)
        dilate(blend_masks[i], blend_masks[i], element);
    // calculate boundary
    for (int iter = 1; iter < stitch_num; iter++) {
        utils::ImageUtils::FindContour(blend_masks[iter], boundaries[iter - 1]);
        unique(boundaries[iter - 1].begin(), boundaries[iter - 1].end());
    }
    // calculate seam
    //cv::Mat tmp_mask(pano_h, pano_w, CV_8U);
    cv::Mat tmp_mask = blend_masks[0].clone();
    const int dx[] = { 0, 1, 1, 1, 0, -1, -1, -1 };
    const int dy[] = { 1, 1, 0, -1, -1, -1, 0, 1 };
    for (int i = 1; i < blend_masks.size(); i++) {
        for (int k = 0; k < boundaries[i - 1].size(); k++) {
            cv::Point pt = boundaries[i - 1][k];
            int ty = pt.y, tx = (pt.x) % pano_w, tx2 = pt.x;
            if (ty < 0 || ty >= pano_h || tx2 < 0 || tx2 >= 2 * pano_w) continue;
            if (blend_masks[i].at<uchar>(ty, tx2) > 0 && (tmp_mask.at<uchar>(ty, tx) > 0)) {
                seams_idx[i - 1].push_back(k);
            }
        }
        for (int row = 0; row < pano_h; row++) {
            for (int col = 0; col < pano_w * 2; col++) {
                uchar cur = blend_masks[i].at<uchar>(row, col);
                if (cur == 0) continue;
                tmp_mask.at<uchar>(row, col % pano_w) = 255;
            }
        }
        /*cv::Mat tmp_show = cv::Mat::zeros(pano_h, pano_w, CV_8U);
        for (int k = 0; k < seams_idx[i - 1].size(); k++) {
            int idx = seams_idx[i - 1].at(k);
            cv::Point p = boundaries[i - 1][idx];
            tmp_show.at<uchar>(p.y, p.x % pano_w) = 255;
        }
        ScaledShow("show", tmp_show);*/
    }
    blend_masks.erase(blend_masks.begin());
}

int main(int argc, char* argv[]) {
    int stitch_num = 6;

    std::vector<cv::Mat> images(stitch_num), pano_masks(stitch_num);
    for (int iter = 0; iter < stitch_num; iter++) {
        std::ostringstream str;
        str << RESOURCES_DIR << iter << ".png";
        images[iter] = cv::imread(str.str(), 1);
        str.str(std::string());
        str << RESOURCES_DIR << "MBB.Mask." << iter << ".png";
        pano_masks[iter] = cv::imread(str.str(), 0);
    }

    int pano_h = pano_masks[0].rows;
    int pano_w = pano_masks[0].cols;
    // calculate boundary
    std::vector<std::vector<cv::Point> > boundaries(stitch_num - 1);
    std::vector<std::vector<int> > seams_idx(stitch_num - 1);
    std::vector<cv::Mat> blend_masks;
    CalculateBoundary(pano_masks, blend_masks, boundaries, seams_idx);

    MvcBlender blender(blend_masks, boundaries, seams_idx, pano_w, pano_h);

    cv::Mat result;
    blender.Blend(images, result);
    utils::ImageUtils::ScaledShow("result", result);

    return EXIT_SUCCESS;
}
