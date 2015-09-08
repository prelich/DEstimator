# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
#
# Linux fixup/packaging script for "make install"

include(/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/../../cmake/BundleUtils-Linux.cmake)
if(FIXUP_BINARY)
    set(bin "test_destimator")
elseif()
    set(bin "")
endif()
set(lib "libdestimator.so")

fixup_linux_libs("${bin}" "${lib}")
