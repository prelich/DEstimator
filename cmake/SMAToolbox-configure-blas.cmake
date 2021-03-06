# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
# Include this file to configure BLAS integration with armadillo

#It appears this is the option required for use with matlab's internal lapack/blas
#which is through the intel MKL uses the long interface
add_definitions(-DARMA_BLAS_LONG)
if(WIN32)
    #These compiled libraries are provided by armadillo
    find_library(BLAS_LIBRARIES blas_win64_MT.dll)
    set(W64_DLLS ${W64_DLLS} blas_win64_MT.dll)
else()
    find_package(BLAS)
endif()
