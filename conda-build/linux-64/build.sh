#!/usr/bin/env bash
echo $PREFIX
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j8 ims-bin && cp ims $PREFIX/bin/ims
