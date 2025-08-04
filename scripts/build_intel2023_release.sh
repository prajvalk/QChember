#!/bin/bash

#SBATCH --job-name=newscf-intel23-release
#SBATCH --times=01:00:00
#SBATCH --output=newscf-intel23-release-%j.out

set -e

. /opt/intel/oneapi/2023/setvars.sh

mkdir build
pushd build
cmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc -DCMAKE_BUILD_TYPE=Release
make -j
make -j test
popd