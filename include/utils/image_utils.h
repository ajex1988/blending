#ifndef INCLUDE_UTILS_IMAGEUTILS_H_
#define INCLUDE_UTILS_IMAGEUTILS_H_

#include <string>
#include <opencv2\opencv.hpp>
namespace utils {

class ImageUtils {
public:
	static void ScaledShow(std::string, const cv::Mat&, double = 1);
	static void FindContour(const cv::Mat&, std::vector<cv::Point>&, int = CV_CHAIN_APPROX_NONE);
};

}

#endif // INCLUDE_UTILS_IMAGEUTILS_H_
