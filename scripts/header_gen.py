#!/usr/bin/env python

#
# Copyright (c) 2018 loyave
#

import filecmp
import inspect
import shutil
import os
import utils

dirname = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))


def main():
    include_dir = os.path.join(dirname, "../include/blending")
    filenames = utils.get_all_files(include_dir, ["*.h"])
    filenames.sort()
    header = os.path.join(dirname, "../include/blending/blending.h")
    header_tmp = header + ".tmp"
    with open(header_tmp, "w") as header_file:
        header_file.write("""// Copyright (c) 2018 loyave
""")
        header_file.write("#ifndef INCLUDE_BLENDING_BLENDING_H_\n")
        header_file.write("#define INCLUDE_BLENDING_BLENDING_H_\n")
        for filename in filenames:
            if not filename.endswith("-inl.h"):
                line = "#include <blending/%s>\n" % os.path.basename(filename)
                header_file.write(line)
        header_file.write("#endif  // INCLUDE_BLENDING_BLENDING_H_\n")
    if not filecmp.cmp(header, header_tmp):
        shutil.move(header_tmp, header)
    else:
        os.remove(header_tmp)


if __name__ == "__main__":
    main()
