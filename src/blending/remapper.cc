#include <blending/remapper.h>
#include <blending/transform.h>

#include <algorithm>

using namespace blending;

Remapper::Remapper(std::string script_name) : ps(script_name) {
    stitch_num = ps.num_images;
    panorama_roi.width = ps.panorama.w;
    panorama_roi.height = ps.panorama.h;
    // initial transforms
    initial_transforms();
    // initial remap table
    remap_table_x = new double*[stitch_num];
    remap_table_y = new double*[stitch_num];
    for (size_t idx = 0; idx < stitch_num; idx++) {
        remap_table_x[idx] = new double[panorama_roi.area()];
        remap_table_y[idx] = new double[panorama_roi.area()];
    }
    initial_remap_table();
}

Remapper::~Remapper() {
    for (size_t idx = 0; idx < stitch_num; idx++) {
        delete[] remap_table_x[idx];
        delete[] remap_table_y[idx];
    }
    delete[] remap_table_x;
    delete[] remap_table_y;
}

std::vector<cv::Mat> Remapper::get_masks() {
    return masks;
}

// Cubic polynomial with parameter A
// A = -1: sharpen; A = - 0.5 homogeneous
// make sure x >= 0
#define	A	(-0.75)
// 0 <= x < 1
inline double Remapper::Cubic01(double x) {
    return	((A + 2.0)*x - (A + 3.0))*x*x + 1.0;
}
// 1 <= x < 2
inline double Remapper::Cubic12(double x) {
    return	((A * x - 5.0 * A) * x + 8.0 * A) * x - 4.0 * A;
}
#undef A

#define 	NNEIGHBOR(x, a , NDIM )											\
	a[0] = 1.0;

#define 	BILINEAR(x, a, NDIM )											\
	a[1] = x;														\
	a[0] = 1.0 - x;

// Unused; has been replaced by 'CUBIC'.
#define 	POLY3( x, a , NDIM )											\
	a[3] = (x * x - 1.0) * x / 6.0;								\
	a[2] = ((1.0 - x) * x / 2.0 + 1.0) * x; 						\
	a[1] = ((1.0 / 2.0  * x - 1.0) * x - 1.0 / 2.0) * x + 1.0;		\
	a[0] = ((-1.0 / 6.0 * x + 1.0 / 2.0) * x - 1.0 / 3.0) * x;

#define		SPLINE16( x, a, NDIM )											\
	a[3] = ((1.0 / 3.0  * x - 1.0 / 5.0) * x - 2.0 / 15.0) * x;		\
	a[2] = ((6.0 / 5.0 - x) * x + 4.0 / 5.0) * x;				\
	a[1] = ((x - 9.0 / 5.0) * x - 1.0 / 5.0) * x + 1.0;		\
	a[0] = ((-1.0 / 3.0 * x + 4.0 / 5.0) * x - 7.0 / 15.0) * x;


#define		CUBIC( x, a, NDIM )									\
	a[3] = Cubic12(2.0 - x);									\
	a[2] = Cubic01(1.0 - x);									\
	a[1] = Cubic01(x);											\
	a[0] = Cubic12(x + 1.0);									\

void Remapper::Remap(const std::vector<cv::Mat>& frames, std::vector<cv::Mat>& remap_frames) {
    remap_frames.resize(stitch_num);
    for (int i = 0; i < stitch_num; i++) {
        remap_frames[i] = cv::Mat::zeros(panorama_roi.size(), CV_8UC3);
    }
    int wd = panorama_roi.width, hd = panorama_roi.height;
    for (int idx = 0; idx < stitch_num; idx++) { // idx : image iterator
        int ws = frames[idx].cols, hs = frames[idx].rows;
        for (int y = 0; y < hd; y++) {
            for (int x = 0; x < wd; x++) {
                double xx = remap_table_x[idx][y * wd + x], yy = remap_table_y[idx][y * wd + x];
                if (xx < 0 || yy < 0) {
                    xx = abs(xx), yy = abs(yy);
                    int xc = std::max(0, std::min((int)floor(xx), ws - 1)), yc = std::max(0, std::min((int)floor(yy), hs - 1));
                    remap_frames[idx].at<cv::Vec3b>(y, x) = frames[idx].at<cv::Vec3b>(yc, xc);
                    continue;
                }
                // interpolation
                int xc = (int)floor(xx), yc = (int)floor(yy);
                double dx = xx - floor(xx), dy = yy - floor(yy);
                int n = 4;
                int n2 = n / 2;
                int ys = yc + 1 - n2; // smallest y-index used for interpolation
                int xs = xc + 1 - n2; // smallest x-index used for interpolation

                double xw[4], yw[4];
                //POLY3(dx, xw, 4);
                //POLY3(dy, yw, 4);
                //SPLINE16(dx, xw, 4);
                //SPLINE16(dy, yw, 4);
                CUBIC(dx, xw, 4);
                CUBIC(dy, yw, 4);

                double weight = 0.0;
                cv::Vec3d sum(0, 0, 0); // b g r
                for (int i = 0; i < n; i++) {
                    int srcy = ys + i;
                    if (srcy < 0 || srcy >= hs) continue;

                    cv::Vec3d sumX(0, 0, 0);
                    for (int j = 0; j < n; j++) {
                        int srcx = xs + j;
                        if (srcx < 0 || srcx >= ws) continue;

                        cv::Vec3b color = frames[idx].at<cv::Vec3b>(srcy, srcx);
                        for (int k = 0; k < 3; k++) sumX[k] += (double)color[k] * xw[j];
                    }
                    for (int k = 0; k < 3; k++) {
                        sum[k] += sumX[k] * yw[i];
                    }
                }
                for (int k = 0; k < 3; k++) {
                    sum[k] = std::max(0.0, std::min((double)UCHAR_MAX, sum[k]));
                    remap_frames[idx].at<cv::Vec3b>(y, x)[k] = (uchar)sum[k];
                }
            }
        }
    }
}

void Remapper::initial_remap_table() {
    std::cout << "Set the frame remapping table : ";
    clock_t start = clock();

    masks.resize(stitch_num);
    for (int i = 0; i < stitch_num; i++) {
        masks[i] = cv::Mat::zeros(panorama_roi.size(), CV_8U);
    }

    int w = ps.panorama.w, h = ps.panorama.h;
    double w2 = (double)w / 2.0 - 0.5, h2 = (double)h / 2.0 - 0.5;
    double dx, dy;
    for (int idx = 0; idx < stitch_num; idx++) {
        int iw = ps.images[idx].w, ih = ps.images[idx].h;
        double sw2 = (double)iw / 2.0 - 0.5, sh2 = (double)ih / 2.0 - 0.5;
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                double y_d = (double)y - h2;
                double x_d = (double)x - w2;
                Transform(idx, x_d, y_d, &dx, &dy);
                dx += sw2;
                dy += sh2;
                if ((dx < (double)iw) && (dx >= 0) && (dy < (double)ih) && (dy >= 0)) {
                    remap_table_x[idx][y * w + x] = dx;
                    remap_table_y[idx][y * w + x] = dy;
                    // set mask
                    if (y < h && x < w)
                        masks[idx].at<uchar>(y, x) = 255;
                }
                else {
                    dx = abs(dx), dy = abs(dy);
                    int xs = (int)floor(dx), ys = (int)floor(dy);
                    dx -= (double)xs;
                    dy -= (double)ys;

                    int xc = xs % (2 * iw), yc = ys % (2 * ih);
                    if (xc >= iw) {
                        xc = 2 * iw - 1 - xc;
                        dx = -dx;
                    }
                    if (yc >= ih) {
                        yc = 2 * ih - 1 - yc;
                        dy = -dy;
                    }
                    // set xc, yc < 0;
                    remap_table_x[idx][y * w + x] = (double)-(dx + (double)xc);
                    remap_table_y[idx][y * w + x] = (double)-(dy + (double)yc);
                }
            }
        }
    }
    std::cout << "Done. " << clock() - start << " ms" << std::endl;
}

void Remapper::initial_transforms() {
    transforms.resize(stitch_num);
    for (size_t idx = 0; idx < stitch_num; idx++) {
        const PanoParam& panorama = ps.panorama;
        const ImageParam& image = ps.images[idx]; // 设置当前图像
        TransformParam& transform = transforms[idx];

        int image_selection_width = image.w;
        int image_selection_height = image.h;

        // 全景图和某路帧对应的视角
        double a = DEG_TO_RAD(image.hfov);
        double b = DEG_TO_RAD(panorama.hfov);

        // 设置变换矩阵
        SetMatrix(-DEG_TO_RAD(image.pitch),
            0.0,
            -DEG_TO_RAD(image.roll),
            transform.mt,
            0);

        transform.distance = ((double)panorama.w) / b;
        // 根据图像类型，计算缩放参数
        switch (image.f) {
        case _rectilinear:
            // calculate distance for this projection
            transform.scale[0] = (double)image_selection_width / (2.0 * tan(a / 2.0)) / transform.distance;
            break;
        case _equirectangular:
        case _panorama:
        case _fisheye_ff:
        case _fisheye_circ:
            transform.scale[0] = ((double)image_selection_width) / a / transform.distance;
            break;
        default:
            std::cerr << "SetMakeParams: Unsupported input image projection" << std::endl;
            break;
        }
        transform.scale[1] = transform.scale[0];
        transform.rot[0] = transform.distance * PI;                                // 180 in screenpoints
        transform.rot[1] = -image.yaw * transform.distance * PI / 180.0;            // rotation angle in screenpoints
        transform.perspect[0] = (void*)(transform.mt);
        transform.perspect[1] = (void*)&(transform.distance);
        // 镜像畸变校正参数
        for (int i = 0; i < 4; i++)
            transform.rad[i] = image.rad[i];
        transform.rad[5] = image.rad[4];
        transform.rad[4] = ((double)image.h) / 2.0;
        // 调节水平竖直偏移
        if (image.horizontal) transform.horizontal = image.hori;
        else transform.horizontal = 0.0;
        if (image.vertical) transform.vertical = image.vert;
        else transform.vertical = 0.0;
    }
}

void Remapper::Transform(int idx, double x, double y, double* dx, double* dy) {
    ImageParam& image = ps.images[idx];
    TransformParam& trans = transforms[idx];

    *dx = x, *dy = y;
    rotate_erect(*dx, *dy, dx, dy, trans.rot); // Rotate equirect. image horizontally
    sphere_tp_erect(*dx, *dy, dx, dy, &(trans.distance)); // Convert spherical image to equirect.
    persp_sphere(*dx, *dy, dx, dy, trans.perspect); // Pers0pective Control spherical Image

    switch (image.f) {
    case _rectilinear:
        rect_sphere_tp(*dx, *dy, dx, dy, &(trans.distance));
        break;
    case _panorama:                                   //  pamoramic image
        pano_sphere_tp(*dx, *dy, dx, dy, &(trans.distance)); // Convert spherical to pano
        break;
    case _equirectangular:                       //  equirectangular image
        erect_sphere_tp(*dx, *dy, dx, dy, &(trans.distance)); // Convert spherical to equirect
        break;
    case _fisheye_circ:
    case _fisheye_ff:
        ; // no conversion needed. It is already in spherical coordinates
        break;
    default:
        std::cerr << "Invalid input projection " << image.f << ". Assumed fisheye." << std::endl;
        exit(0);
    }
    resize(*dx, *dy, dx, dy, trans.scale);
    radial(*dx, *dy, dx, dy, trans.rad);

    if (trans.horizontal != 0.0) {
        horiz(*dx, *dy, dx, dy, &(trans.horizontal));
    }
    if (trans.vertical != 0.0) {
        vert(*dx, *dy, dx, dy, &(trans.vertical));
    }
}
