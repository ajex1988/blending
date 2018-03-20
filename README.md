# Panorama Stitch Tool - Blending

A panorama stitcher based on PTGUI template.

Use PTGUI tool to generate template for static scenes.
Videos are captured by a rig consisted of six gopros.
Different blending algorithms are applied to the videos including feather blending, multiband blending and mean-value coordinate method.

## Requirements
* CGAL >= 4.10
* OpenCV >= 3.2 WITH CUDA 8.0

## Quick Start

You will need CMake to build the code. If you're using Windows, you need Visual Studio 2015 or 2017 in addition to CMake.

First, clone the code:

```
git clone https://github.com/ajex1988/blending.git
cd blending
```

### C++ API

For macOS or Linux:

```
mkdir build && cd build && cmake .. && make
```

For Windows:

```
mkdir build
cd build
cmake .. -G"Visual Studio 14 2015 Win64"
MSBuild blending.sln /p:Configuration=Release
```

Now run some examples, such as:

```
bin/mvc_blender_demo
```

## License
