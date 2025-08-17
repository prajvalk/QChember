#!/usr/bin/env bash
set -euo pipefail

echo "=== Build Config ==="
echo "COMPILER=$COMPILER"
echo "BLAS_VENDOR=$BLAS_VENDOR"
echo "OPENMP=$OPENMP"
echo "BUILD_TYPE=$BUILD_TYPE"
echo "===================="

# Update packages
sudo apt-get update -y
sudo apt-get install -y git cmake make wget

# Select compiler
case "$COMPILER" in
  gcc)
    sudo apt-get install -y gcc g++ gfortran
    export CC=gcc
    export CXX=g++
    ;;
  clang)
    sudo apt-get install -y clang lld gfortran
    export CC=clang
    export CXX=clang++
    ;;
  intel)
    wget https://registrationcenter-download.intel.com/akdlm/irc_nas/tec/oneapi-latest/linux/oneapi-latest.tar.gz
    tar -xzf oneapi-latest.tar.gz
    ./oneapi*/bootstrapper --action install --components intel-icc,ifort,mkl --silent
    source /opt/intel/oneapi/setvars.sh
    export CC=icc
    export CXX=icpc
    ;;
  amd)
    wget https://download.amd.com/developer/eula/aocc/aocc-5-0/aocc-compiler-5.0.0.tar
    tar -xf aocc-compiler.tar
    ./aocc*/install.sh --accept
    source /opt/AMD/aocc*/setenv.sh
    export CC=clang
    export CXX=clang++
    ;;
esac

# Install BLAS/LAPACK
case "$BLAS_VENDOR" in
  openblas)
    sudo apt-get install -y libopenblas-dev liblapack-dev
    ;;
  mkl)
    # Intel MKL installed already with oneAPI above
    ;;
  blis)
    sudo apt-get install -y libblis-dev libflame-dev
    ;;
esac

# Configure OpenMP
if [[ "$OPENMP" == "ON" ]]; then
  OMP_FLAG="-DUSE_OPENMP=ON"
else
  OMP_FLAG="-DUSE_OPENMP=OFF"
fi

# Build
mkdir -p build
cmake -B build \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  $OMP_FLAG
cmake --build build -j
