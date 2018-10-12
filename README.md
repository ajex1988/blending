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

```bash
mkdir build && cd build && cmake .. && make
```

For Windows:

```sh
mkdir build
cd build
cmake .. -G"Visual Studio 14 2015 Win64"
MSBuild blending.sln /p:Configuration=Release
```

Now run some examples, such as:

```sh
bin/mvc_blender_demo
```


## Resources
Click [here](http://cg.cs.tsinghua.edu.cn/blending/files/resources.zip) to download the resources file.
It contains orignal images to be stitched, a PTStitcher script file generate by PTGui tool, and the mask images generate by seam finder demo program. 

You should unzip the file to */resources/* folder to run the demo programs under */src/examples*.
