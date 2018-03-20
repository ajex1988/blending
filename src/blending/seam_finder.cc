#include <blending/seam_finder.h>
#include <utils/image_utils.h>>

#include <ctime>
#include <opencv2/stitching/detail/seam_finders.hpp>

using namespace blending;

SeamFinder::SeamFinder() {
}


SeamFinder::~SeamFinder() {
}

void SeamFinder::find(std::vector<cv::Mat>& masks, int _w, int _h, bool expand) {
	imNum = masks.size();
	w = _w, h = _h;
	for (int i = 0; i < imNum - 1; i++) {
		for (int j = i + 1; j < imNum; j++) {
			clock_t s = clock(), e;
			findInPair(masks[i], masks[j], expand);
			e = clock();
			std::cout << "  generate masks for " << i << ' ' << j << ' ' << (e - s) / 1000.0 << " s." << std::endl;
		}
	}
	return;
}

void SeamFinder::findInPair(cv::Mat m1, cv::Mat m2, bool expand) {
	cv::Mat fm1, fm2;
	fm1.create(h, 2 * w, CV_8U);
	fm2.create(h, 2 * w, CV_8U);
	cv::Range r1 = getRange(m1), r2 = getRange(m2);
	cv::Range r;
	int xmin = MAX(r1.start, r2.start), xmax = MIN(r1.end, r2.end);
	if (xmin <= xmax) {
		if ((r1.end >= w && (r1.end != 2 * w - 1)) || (r2.end >= w && (r2.end != 2 * w - 1))) {
			r.start = w / 2;
			r.end = w / 2 + w - 1;
		}
		else {
			r.start = 0;
			r.end = w - 1;
		}
	}
	else {
		if (r1.end < w && r2.end < w) return;
		if (r1.end < w) r1 = r1 + w - 1;
		if (r2.end < w) r2 = r2 + w - 1;
		xmin = MAX(r1.start, r2.start), xmax = MIN(r1.end, r2.end);
		if (xmin <= xmax) {
			r.start = w / 2;
			r.end = w / 2 + w - 1;
		}
		else {
			return;
		}
	}
	middleCopy(fm1, m1, r, COPY2FULL);
	middleCopy(fm2, m2, r, COPY2FULL);
	middleCut(fm1, fm2, r, expand);
	middleCopy(fm1, m1, r, COPY2SMALL);
	middleCopy(fm2, m2, r, COPY2SMALL);
}

cv::Range SeamFinder::getRange(cv::Mat m) {
	w = m.cols, h = m.rows;

	int cl = 0, cr = 0, ct = 0, cb = 0;
	int xmin = m.cols * 2, xmax = -1;
	for (int i = 0; i < m.rows; i++) {
		if (m.at<uchar>(i, 0) > 0) cl++;
		if (m.at<uchar>(i, m.cols - 1) > 0) cr++;
	}
	for (int i = 0; i < m.cols; i++) {
		if (m.at<uchar>(0, i) > 0) ct++;
		if (m.at<uchar>(m.rows - 1, 0) > 0) cb++;
	}
	if ((ct == w) || cb == w) { // full-ceiling
		xmin = 0;
		xmax = 2 * w - 1;
	}
	else if (cl > 0 && cr > 0) {
		for (int i = 0; i < m.rows; i++) {
			for (int j = 0; j < m.cols; j++) {
				if (m.at<uchar>(i, j) == 0) continue;
				int x = (j < m.cols / 2) ? j + m.cols : j;
				if (x > xmax) xmax = x;
				if (x < xmin) xmin = x;
			}
		}
	}
	else {
		for (int i = 0; i < m.rows; i++) {
			for (int j = 0; j < m.cols; j++) {
				if (m.at<uchar>(i, j) == 0) continue;
				if (j > xmax) xmax = j;
				if (j < xmin) xmin = j;
			}
		}
	}
	return cv::Range(xmin, xmax);
}

void SeamFinder::middleCut(cv::Mat m1, cv::Mat m2, cv::Range r, bool expand) {
	int xmin = r.start, xmax = r.end;
	std::vector<cv::Point> c, ta, tb;
	for (int i = 0; i < m1.rows; i++) {
		for (int j = xmin; j <= xmax; j++) {
			uchar p1 = m1.at<uchar>(i, j), p2 = m2.at<uchar>(i, j);
			if (p1 > 0 && p2 > 0) {
				c.push_back(cv::Point(j, i));
			}
			else if (p1 > 0) {
				ta.push_back(cv::Point(j, i));
			}
			else if (p2 > 0) {
				tb.push_back(cv::Point(j, i));
			}
		}
	}
	for (int i = 0; i < c.size(); i++) {
		m1.at<uchar>(c[i].y, c[i].x) = 0;
		m2.at<uchar>(c[i].y, c[i].x) = 0;
	}
	std::vector<cv::Point> a, b;
	utils::ImageUtils::FindContour(m1, a);
	utils::ImageUtils::FindContour(m2, b);
	for (int i = 0; i < c.size(); i++) {
		double suma = 0, sumb = 0;
		for (int j = 0; j < a.size(); j++) {
			cv::Point d = a[j] - c[i];
			suma += 1 / sqrt(d.ddot(d));
		}
		for (int j = 0; j < b.size(); j++) {
			cv::Point d = b[j] - c[i];
			sumb += 1 / sqrt(d.ddot(d));
		}
		double eps = (suma + sumb) / 50;
		if (!expand) {
			if (suma  >= sumb) {
				m1.at<uchar>(c[i].y, c[i].x) = 255;
			}
			else {
				m2.at<uchar>(c[i].y, c[i].x) = 255;
			}
		}
		else {
			if (suma + eps >= sumb) {
				m1.at<uchar>(c[i].y, c[i].x) = 255;
			}
			if (sumb + eps >= suma){
				m2.at<uchar>(c[i].y, c[i].x) = 255;
			}
		}
	}
}

void SeamFinder::middleCopy(cv::Mat fmask, cv::Mat mask, cv::Range r, int type) {
	if (type == COPY2SMALL) {
		if (r.end < w) {
			fmask(cv::Rect(0, 0, w, h)).copyTo(mask);
		}
		else if (r.start >= w) {
			fmask(cv::Rect(w, 0, w, h)).copyTo(mask);
		}
		else {
			fmask(cv::Rect(w, 0, w / 2, h)).copyTo(mask(cv::Rect(0, 0, w / 2, h)));
			fmask(cv::Rect(w / 2, 0, w / 2, h)).copyTo(mask(cv::Rect(w / 2, 0, w / 2, h)));
		}
	}
	else if (type == COPY2FULL) {
		if (r.end < w)
			mask.copyTo(fmask(cv::Rect(0, 0, w, h)));
		else {
			mask(cv::Rect(0, 0, w / 2, h)).copyTo(fmask(cv::Rect(w, 0, w / 2, h)));
			mask(cv::Rect(w / 2, 0, w / 2, h)).copyTo(fmask(cv::Rect(w / 2, 0, w / 2, h)));
		}
	}
}

void SeamFinder::findVoronoi(std::vector<cv::Mat>& masks, int w, int h) {
	std::vector<cv::Size> sizes;
	std::vector<cv::Point> corners;
    std::vector<cv::UMat> u_masks;
	for (int i = 0; i < masks.size(); i++) {
		sizes.push_back(cv::Size(w, h));
		corners.push_back(cv::Point(0, 0));
        u_masks.push_back(masks[i].getUMat(cv::ACCESS_RW));
	}
    cv::detail::VoronoiSeamFinder vsf;
	vsf.find(sizes, corners, u_masks);
	int dilation_size = 10;
	cv::Mat element = getStructuringElement(cv::MORPH_ELLIPSE,
		cv::Size(2 * dilation_size + 1, 2 * dilation_size + 1),
		cv::Point(dilation_size, dilation_size));
	for (int i = 0; i < masks.size(); i++) {

		//cv::Mat tmp;
		//masks[i].copyTo(tmp);
		dilate(masks[i], masks[i], element);
		//compare(tmp, masks[i], tmp, CMP_NE);
		//PanoStitch::showImage(tmp, 0.6);
		/*GaussianBlur(masks[i], masks[i], cv::Size(3, 3), 0, 0);
		compare(masks[i], 0, masks[i], CMP_NE);*/
	}
}
