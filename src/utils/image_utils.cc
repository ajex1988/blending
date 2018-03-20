#include <utils/image_utils.h>

using namespace utils;

/// show scaled images;
void ImageUtils::ScaledShow(std::string title, const cv::Mat& img, double scale) {
    cv::Mat tmp;
	img.copyTo(tmp);
	cv::resize(tmp, tmp, cv::Size(), scale, scale);
	cv::imshow(title, tmp);
	cv::waitKey();
}
/// find contours. Fix Opencv::findContours image-border clipping bug.
/// avoid border clipping in OpenCV findContours
/// choose the largest contour
void ImageUtils::FindContour(const cv::Mat & img, std::vector<cv::Point>& contours, int approx) {
	contours.clear();
	cv::Mat tmp = cv::Mat::zeros(img.rows + 2, img.cols + 2, img.type());
	img.copyTo(tmp(cv::Rect(1, 1, img.cols, img.rows)));
	std::vector<std::vector<cv::Point> > _contours;
	std::vector<cv::Vec4i> hierarchy;
	findContours(tmp, _contours, hierarchy, CV_RETR_EXTERNAL, approx, cv::Point(-1, -1));
	int idx = 0, max = 0;
	for (int i = 0; i < _contours.size(); i++) {
		if (_contours[i].size() > max) {
			max = _contours[i].size();
			idx = i;
		}
	}
	contours = _contours[idx];
}
