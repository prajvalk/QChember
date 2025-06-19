SCRIPT_VERSION = "0.1.0"
PROJECT_VERSION = ""

DEFAULT_FLAGS="-fpermissive"
DEBUG_FLAGS=" -g3" # Modify Later
RELEASE_FLAGS=" -O3"
CXX_LIB_FLAGS=" -fPIC -shared"

AVX2_FLAG=" -mavx2"
AVX512_FLAGS=" -mavx512f -mavx512cd -mavx512dq -mavx512bw -mavx512vl"

with open("VERSION", "r") as file:
    PROJECT_VERSION = file.readline()

output_data = ""

def get_index(options, cmp):
    i = 0
    for opt in options:
        if opt == cmp:
            return i
        else:
            i=i+1
    return -1

def get_user_input(prompt, options, default_index):
    while True:
        inp = input(prompt+" "+str(options)+" ["+options[default_index]+"]: ")
        if inp == "" or inp == None:
            return default_index
        if inp in options:
            return get_index(options, inp)
        else:
            print("Available choices:", options)

def get_compiler_flags(compiler, buildtype, simd_level, openmp):
    flags=DEFAULT_FLAGS
    if buildtype == "release":
        flags += RELEASE_FLAGS
        if compiler == "g++" or compiler == "clang++":
            flags += " -flto"
        else:
            flags += " -ipo"
    else:
        flags += DEBUG_FLAGS
        if compiler == "g++" or compiler == "clang++":
            flags += " -Og"
        else:
            flags += " -O0"
    if simd_level == "native":
        if compiler == "g++" or compiler == "clang++":
            flags += " -march=native"
        else:
            flags += " -xHost"
    elif simd_level == "avx2":
        flags += AVX2_FLAG
    elif simd_level == "avx512":
        flags += AVX512_FLAGS
    else:
        if compiler == "g++":
            flags += " -mno-sse -mno-avx -fno-tree-vectorize"
        elif compiler == "clang++":
            flags += " -mno-sse -mno-avx -fno-vectorize -fno-slp-vectorize"
        else:
            flags += " -no-vec -mno-sse -mno-avx"
    if openmp == "openmp":
        flags += " -fopenmp"
    else: 
        flags += "-fno-openmp"
    return flags
    

build_flags=["release", "debug"]
simd_flags=["native", "avx2", "avx512", "none"]
openmp_flags=["yes", "no"]
compilers=["g++", "clang++", "icpx"]
blas_type=["none", "netlib", "openblas", "flexiblas", "intel_mkl", "amd_blis"]

print("new_scf Python Build Bootstrapping System (pybb)")
print("Version", SCRIPT_VERSION)
print("Mozilla Public License 2.0, see LICENSE.md for the full license terms.")
print("Copyright (C) 2025, Prajval K")
print()

user_build=get_user_input("Build Type:", build_flags, 0)
user_simd=get_user_input("CPU Intrinsics Level:", simd_flags, 0)
user_openmp=get_user_input("Use OpenMP for parallelization?", openmp_flags, 0)
user_compiler=get_user_input("C++ Compiler:", compilers, 0)
user_blas=get_user_input("BLAS Type:", blas_type, 0)

user_build=build_flags[user_build]
user_simd=simd_flags[user_simd]
if openmp_flags[user_openmp] == "yes":
    user_openmp="openmp"
else:
    user_openmp="serial"
user_compiler=compilers[user_compiler]
user_blas=blas_type[user_blas]

config_tag=f"newscf-{PROJECT_VERSION}-{user_build}-{user_simd}-{user_openmp}-{user_blas}"

print()
print("Build Tag:", config_tag)

CXX=user_compiler
CXX_FLAGS=get_compiler_flags(CXX, user_build, user_simd, user_openmp)

print(CXX_FLAGS)