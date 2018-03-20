#include <blending/blending.h>

#define APP_NAME "seam_finder_demo"

#include <direct.h>

using namespace blending;

int main(int argc, char* argv[]) {
    const std::string script_file = RESOURCES_DIR "script.txt";
    const std::string mask_output_dir = APP_NAME "_masks/";
    _mkdir(mask_output_dir.c_str());

    std::vector<cv::Mat> remap_masks;
    Remapper pr(script_file);
    remap_masks = pr.get_masks();

    int stitch_num = remap_masks.size();
    int pano_w, pano_h;
    if (stitch_num > 0) {
        pano_w = remap_masks[0].cols;
        pano_h = remap_masks[0].rows;
    }
    else {
        std::cerr << "Error load masks" << std::endl;
        exit(EXIT_FAILURE);
    }

    SeamFinder sf;
    sf.find(remap_masks, pano_w, pano_h, false);

    for (int iter = 0; iter < stitch_num; iter++) {
        std::ostringstream str;
        str << mask_output_dir << iter << ".png";
        cv::imwrite(str.str(), remap_masks[iter]);
    }

    return EXIT_SUCCESS;
}