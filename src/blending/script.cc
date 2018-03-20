#include <blending/script.h>
#include <blending/transform.h>

#include <fstream>
#include <sstream>
#include <iostream>

using namespace blending;

Script::Script(std::string script_name) {
    LoadScript(script_name);
}

Script::~Script() {
}

void Script::LoadScript(const std::string& script_name) {
    std::ifstream fin(script_name);
    if (!fin.is_open()) {
        std::cerr << "Failed to open file : " << script_name << std::endl;
        exit(EXIT_FAILURE);
    }
    num_images = 0;
    while (!fin.eof()) {
        std::string line;
        getline(fin, line);
        switch (line[0]) {
        case 'p':
            ReadPanoLine(line);
            break;
        case 'o':
            images.push_back(ReadImageLine(line));
            num_images++;
            break;
        case 'i': // ignore
        default:
            break;
        }
    }
    fin.close();
}

void Script::ReadPanoLine(const std::string& line) {
    std::istringstream istr(line);
    std::string str;
    istr >> str; // read 'p' charactor
    while (!istr.eof()) {
        istr >> str;
        switch (str[0]) {
        case 'w':
            sscanf(str.c_str(), "w%d", &panorama.w);        // read width
            break;
        case 'h':
            sscanf(str.c_str(), "h%d", &panorama.h);        // read height
            break;
        case 'v':
            sscanf(str.c_str(), "v%lf", &panorama.hfov);    // read fov
            break;
        default:
            break;
        }
    }
}

// read 'o' lines;
// ignore '+buff', '-buff', 'u'.etc(obsolete) and 'o'
ImageParam Script::ReadImageLine(const std::string& line) {
    std::istringstream istr(line);
    std::string str;
    istr >> str; // read 'o' charactor
    ImageParam image_param;
    image_param.panorama = &(this->panorama);
    while (!istr.eof()) {
        istr >> str;
        switch (str[0]) {
        case 'w':
            sscanf(str.c_str(), "w%d", &(image_param.w));
            break;
        case 'h':
            sscanf(str.c_str(), "h%d", &(image_param.h));
            break;
        case 'f':
            sscanf(str.c_str(), "f%d", &(image_param.f));
            break;
        case 'v':
            sscanf(str.c_str(), "v%lf", &(image_param.hfov));
            break;
        case 'r':
            sscanf(str.c_str(), "r%lf", &(image_param.roll));
            break;
        case 'p':
            sscanf(str.c_str(), "p%lf", &(image_param.pitch));
            break;
        case 'y':
            sscanf(str.c_str(), "y%lf", &(image_param.yaw));
            break;
        case 'a':
            sscanf(str.c_str(), "a%lf", &(image_param.rad[3]));
            image_param.radial = true;
            break;
        case 'b':
            sscanf(str.c_str(), "b%lf", &(image_param.rad[2]));
            image_param.radial = true;
            break;
        case 'c':
            sscanf(str.c_str(), "c%lf", &(image_param.rad[1]));
            image_param.radial = true;
            break;
        case 'd':
            sscanf(str.c_str(), "d%lf", &(image_param.hori));
            image_param.horizontal = true;
            break;
        case 'e':
            sscanf(str.c_str(), "e%lf", &(image_param.vert));
            image_param.vertical = true;
            break;
        case 'S':
            sscanf(str.c_str(), "S%d,%d,%d,%d", &(image_param.S[0]), &(image_param.S[1]), &(image_param.S[2]), &(image_param.S[3]));
            image_param.select = true;
            break;
        case 'C':
            sscanf(str.c_str(), "C%d,%d,%d,%d", &(image_param.C[0]), &(image_param.C[1]), &(image_param.C[2]), &(image_param.C[3]));
            image_param.crop = true;
            break;
        default:
            break;
        }
    }
    image_param.rad[0] = 1.0 - (image_param.rad[1] + image_param.rad[2] + image_param.rad[3]);
    SetCorrectionRadius(image_param.rad);
    return image_param;
}
