#ifndef INCLUDE_BLENDING_SCRIPT_H_
#define INCLUDE_BLENDING_SCRIPT_H_

#include <blending/common.h>

#include <string>
#include <vector>

namespace blending {

class Script {
public:
    Script(std::string);
    ~Script();

    void LoadScript(const std::string&);
private:
    void ReadPanoLine(const std::string&);
    ImageParam ReadImageLine(const std::string&);
public:
    int num_images;
    PanoParam panorama;
    std::vector<ImageParam> images;
};

}

#endif // INCLUDE_BLENDING_SCRIPT_H_
