#!/bin/bash

set -e
set -x

CXX=g++
CXX_FLAGS="-O3 -flto -g3"
CXX_LIB_FLAGS="-fPIC -shared"
PROJECT_ROOT=$(pwd)

rm -rf ./build
mkdir -p build
mkdir -p build/lib
mkdir -p build/tests
BUILD_DIR=$(pwd)/build
LIB_DIR=$BUILD_DIR/lib
TEST_DIR=$BUILD_DIR/tests

INCLUDE_FLAGS=""

##################
# new_scf/common #
##################
INCLUDE_FLAGS="-I$PROJECT_ROOT/common/include $INCLUDE_FLAGS"

####################
# new_scf/blasbind #
####################
pushd blasbind
BLASBIND_INCLUDE="-I$(pwd)/include"
BLASBIND_LINK_LIB=""
BLASBIND_OPT_FLAGS="-fopenmp -march=native"

$CXX $CXX_FLAGS $CXX_LIB_FLAGS $INCLUDE_FLAGS $BLASBIND_INCLUDE $BLASBIND_LINK_LIB $BLASBIND_OPT_FLAGS src/*.cpp -o $LIB_DIR/libblasbind.so

# Tests
$CXX $CXX_FLAGS $INCLUDE_FLAGS $BLASBIND_INCLUDE $BLASBIND_LINK_LIB $BLASBIND_OPT_FLAGS tests/gemm_test.cpp -o $TEST_DIR/blasbind_gemm_test -L$LIB_DIR -lblasbind
LD_LIBRARY_PATH=$LIB_DIR $TEST_DIR/blasbind_gemm_test

popd

#####################
# new_scf/libmatrix #
#####################
pushd libmatrix

LIBMATRIX_INCLUDE="-I$(pwd)/include"
LIBMATRIX_LINK_LIB=""
LIBMATRIX_OPT_FLAGS="-fopenmp -march=native"
LIBMATRIX_SRC="src/*.cpp"

# Compile Lib
$CXX $CXX_FLAGS $CXX_LIB_FLAGS $INCLUDE_FLAGS $LIBMATRIX_INCLUDE $LIBMATRIX_LINK_LIB $LIBMATRIX_OPT_FLAGS $LIBMATRIX_SRC -o $LIB_DIR/libmatrix.so

# Compile Tests
$CXX $CXX_FLAGS $INCLUDE_FLAGS $LIBMATRIX_INCLUDE $LIBMATRIX_LINK_LIB $LIBMATRIX_OPT_FLAGS tests/matrix_test.cpp -o $TEST_DIR/matrix_test -L$LIB_DIR -lmatrix

# Run Tests
LD_LIBRARY_PATH=$LIB_DIR $TEST_DIR/matrix_test

popd