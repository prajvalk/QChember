#!/bin/bash

#SBATCH --job-name=newscf-gcc-release
#SBATCH --times=01:00:00
#SBATCH --output=newscf-gcc-release-%j.out

set -e

mkdir build
pushd build
cmake .. -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release
make -j
make -j test
popd