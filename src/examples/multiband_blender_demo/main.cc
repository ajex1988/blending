#include <blending/blending.h>
#include <utils/image_utils.h>

#define APP_NAME "multiband_blender_demo"

using namespace blending;

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

    MultiBandBlender blender(pano_masks);

    cv::Mat result;
    blender.Blend(images, result);
    utils::ImageUtils::ScaledShow("result", result);

    return EXIT_SUCCESS;
}
