# new_scf : An optimized* quantum chemisty calculation library
*WIP

## pybb (bootstrap.py) -- Python Build Bootstrap Toolkit
Creates a ```checkhost.sh``` script to verify host capabilities for chosen build configuration, along with ```build.sh``` to compile the library, and a ```run_tests.sh``` script to run the tests (if selected).

## blasbind
Configurable BLAS API that ships OpenMP (AVX512 and AVX2 capable) accelerated routines by default and is compatible to link with Netlib FBLAS, OpenBLAS, FlexiBLAS API, Intel oneAPI MKL FBLAS, AMD BLIS if detected.

Commons will ship a LAPACKE interface header, which will use Netlib FORTRAN LAPACK, OpenBLAS LAPACK, or MKL or AMD libFLAME depending on the build configuration.

## libmatrix
Lightweight optimized OpenMP Matrix library (AVX512 and AVX2 capable). Templated by default, but will use accelerated routines if standard types are used.

## commons
#### testing_api.hpp
#### logging_api.hpp
#### benchmark_api.hpp
#### lapacke.h

## Planned Capabilities
Sorted according to priority
* (Restricted) Hartree-Fock Calculations for Atoms (atomscf)
* (Unrestricted) Hartree-Fock Calculations for atoms (atomscf)
* (Restricted) Hartree-Fock Calculations for Molecules (molscf)
* MP2 Calculations (mp2)
* CUDA-acceleration for HF and MP2
* Fully incore-GPU computation routines
* __binary128 SIMD intrinsics (AVX512 only) acceleration for libmatrix and blasbind
* GNU MPFR intrinsic support for libmatrix, blasbind

#
Mozilla Public License 2.0 (see LICENSE.md)

Copyright (C) 2025, Prajval K
