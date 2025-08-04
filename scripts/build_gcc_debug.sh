#!/bin/bash

#SBATCH --job-name=newscf-gcc-debug
#SBATCH --times=01:00:00
#SBATCH --output=newscf-gcc-debug-%j.out

set -e

mkdir build
pushd build
cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CMAKE_BUILD_TYPE=Debug
make -j
make -j test
popd