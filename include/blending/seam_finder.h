#ifndef INCLUDE_BLENDING_SEAMFINDER_H_
#define INCLUDE_BLENDING_SEAMFINDER_H_

#include <opencv2\opencv.hpp>

namespace blending {
class SeamFinder {
public:
	SeamFinder();
	~SeamFinder();

	enum {
		COPY2SMALL,
		COPY2FULL
	};

	void find(std::vector<cv::Mat>&, int, int, bool = true);
	void findInPair(cv::Mat, cv::Mat, bool);
	void findVoronoi(std::vector<cv::Mat>&, int, int);
private:
	cv::Range getRange(cv::Mat);
	void middleCut(cv::Mat, cv::Mat, cv::Range, bool);
	void middleCopy(cv::Mat, cv::Mat, cv::Range, int);
	void combine(cv::Mat, cv::Mat, uchar v);

	int imNum, w, h;
};

}

#endif // INCLUDE_BLENDING_SEAMFINDER_H_
