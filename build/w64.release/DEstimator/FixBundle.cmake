# SMAToolbox CMake build system
# Mark J. Olah (mjo@cs.unm.edu)
# 03-2014
#
# Windows fixup/packaging script for "make install"

include(/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/../../cmake/BundleUtils-W64.cmake)
fixup_w64_libs("../../../mex/mex.w64" "libstdc++-6.dll;libgcc_s_seh-1.dll;libgfortran-3.dll;libgomp-1.dll;libwinpthread-1.dll;blas_win64_MT.dll" "/home/mjo/mxe/usr/x86_64-w64-mingw32.shared/lib;/home/mjo/mxe/usr/x86_64-w64-mingw32.shared/bin;/home/mjo/mxe/usr/bin;/home/mjo/mxe/usr/lib;/home/mjo/mxe/usr/lib/gcc/x86_64-w64-mingw32.shared/4.9.2/;/home/mjo/usr/x86_64-w64-mingw32/lib")
