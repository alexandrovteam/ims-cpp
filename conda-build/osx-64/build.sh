#!/usr/bin/env bash

echo $PREFIX
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=gcc-5 -DCMAKE_CXX_COMPILER=g++-5
make -j8 ims-bin && cp ims $PREFIX/bin/ims
