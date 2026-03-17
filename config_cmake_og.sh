# Auto-Generated on Mon Aug 25 13:54:22 2025
# using library build at /Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs
# on machine Sophias-MacBook-Pro.local
# settings: --nthread=2 --build_umfpack=1 --build_arpack=1

# Setup build and install paths
ROOT_PATH=$(pwd)
BUILD_DIR=$ROOT_PATH/build_release
INSTALL_DIR=$ROOT_PATH/install_release

# Create fresh build directory
rm -rf $BUILD_DIR
mkdir $BUILD_DIR && cd $BUILD_DIR

cmake \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX:PATH=$INSTALL_DIR \
  -DOFT_BUILD_TESTS:BOOL=FALSE \
  -DOFT_BUILD_EXAMPLES:BOOL=TRUE \
  -DOFT_BUILD_PYTHON:BOOL=TRUE \
  -DOFT_BUILD_DOCS:BOOL=FALSE \
  -DOFT_USE_OpenMP:BOOL=TRUE \
  -DOFT_PACKAGE_BUILD:BOOL=FALSE \
  -DOFT_PACKAGE_NIGHTLY:BOOL=TRUE \
  -DOFT_COVERAGE:BOOL=FALSE \
  -DOFT_DEBUG_STACK:BOOL=FALSE \
  -DOFT_PROFILING:BOOL=FALSE \
  -DCMAKE_C_COMPILER:FILEPATH=gcc-15 \
  -DCMAKE_CXX_COMPILER:FILEPATH=g++-15 \
  -DCMAKE_Fortran_COMPILER:FILEPATH=gfortran-15 \
  -DCMAKE_Fortran_FLAGS:STRING="-fallow-argument-mismatch" \
  -DOFT_USE_MPI:BOOL=FALSE \
  -DOFT_METIS_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/metis-5_1_0 \
  -DHDF5_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/hdf5-1_14_6 \
  -DBLAS_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/OpenBLAS-0_3_29 \
  -DLAPACK_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/OpenBLAS-0_3_29 \
  -DBLA_VENDOR:STRING=OpenBLAS \
  -DOFT_ARPACK_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/arpack-ng-3_9_1 \
  -DOFT_FoX_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/fox-4_1_2 \
  -DOFT_UMFPACK_ROOT:PATH=/Users/sophiaguizzo/Documents/Research/Tools/OFT/ext_libs/UMFPACK-6_3_5 \
  /Users/sophiaguizzo/Documents/Research/mug2d/OpenFUSIONToolkit/src
