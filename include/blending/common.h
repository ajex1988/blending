#ifndef INCLUDE_BLENDING_COMMON_H_
#define INCLUDE_BLENDING_COMMON_H_

namespace blending {

enum
{                               // Enumerates for Image.format
    _rectilinear = 0,           // (Standand) FOV (rectilinear) =  2 * arctan (frame size/(focal length * 2))
    _panorama = 1,              // Cylindrical
    _fisheye_circ = 2,          // fisheye-equidistance Circular
    _fisheye_ff = 3,            // fisheye-equidistance Full Frame
    _equirectangular = 4,
    _spherical_cp = 5,          // Fisheye-Horizontal is an image shot with the camera held horizontally. The equator is now in the center of the image.
    _spherical_tp = 6,          // Fisheye-vertical is an image shot with the camera held vertically up.  The panorama is extracted from the circumpherence of the image.
    _mirror = 7,                // convex mirror. This is the reflection of a convex, spherical image. The horizontal field of view is calculated using the formula HFov = 2*arcsin(radius of mirror/radius of curvature of mirror)
    _orthographic = 8,          // fisheye-orthographic FOV  (orthogonal fisheye) = 2 * arcsin (frame size/(focal length *2)
    _cubic = 9,
    _stereographic = 10,        // fisheye stereographic FOV (stereographic fisheye) = 4 * arctan (frame size/(focal length * 4))
    _mercator = 11,
    _trans_mercator = 12,
    _trans_panorama = 13,
    _sinusoidal = 14,
    _lambert = 15,
    _lambertazimuthal = 16,
    _albersequalareaconic = 17,
    _millercylindrical = 18,
    _panini = 19,
    _architectural = 20,
    _equisolid = 21,            // fisheye-equisolid  FOV (equisolid fisheye) = 4 * arcsin (frame size/(focal length * 4))
    _equipanini = 22,
    _biplane = 23,
    _triplane = 24,
    _panini_general = 25,
    _thoby = 26,                // generalizes the model found in modern fisheye lenses. It is
                                // parametrizable but it defaults to the Nikkor 10.5 fisheye lens
    _hammer = 27,
};

struct PanoParam final {                  // 记录生成全景图的参数。
    int w, h;                       // 生成全景图的宽度，高度
    double hfov;                    // 生成全景图的视野范围，1~360°
};

struct ImageParam final {                 // 记录每路视频帧的参数
    int w, h, f;                    // 帧宽度，高度，类型
    int S[4], C[4];                 // 剪裁参数， select或crop
    double hfov;                    // 帧的视野范围，1~360°
    double roll, pitch, yaw;        // 三个方向的旋转角度
    double g, t;
    double rad[6], hori, vert;      // 畸变校正，偏移校正
    bool radial;                    // 记录是否需要畸变校正
    bool horizontal, vertical;      // 记录是否发生偏移
    bool select, crop;              // 记录是否发生裁剪
    PanoParam* panorama;            // 对应的全景图参数
};

struct TransformParam final {             // 记录变换的参数
    double distance;
    double rot[2];
    double scale[2];
    double mt[3][3];
    double rad[6];
    double horizontal, vertical;
    void * perspect[2];
};

}
#endif // INCLUDE_BLENDING_COMMON_H_
