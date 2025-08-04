#!/bin/bash

#SBATCH --job-name=newscf-intel25-release
#SBATCH --times=01:00:00
#SBATCH --output=newscf-intel25-release-%j.out

set -e

. /opt/intel/oneapi/2025/setvars.sh

mkdir build
pushd build
cmake .. -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCMAKE_BUILD_TYPE=Release
make -j
make -j test
popd