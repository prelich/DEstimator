# Install script for directory: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "../../..")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Runtime")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64"
         RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib:/opt/MATLAB/R2014b/bin/glnxa64")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64" TYPE SHARED_LIBRARY FILES "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/DEstimator_Iface.mexa64")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64"
         OLD_RPATH "/opt/MATLAB/R2014b/bin/glnxa64:/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src:"
         NEW_RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib:/opt/MATLAB/R2014b/bin/glnxa64")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/mex/mex.glnxa64/DEstimator_Iface.mexa64")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Runtime")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so"
         RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/libdestimator.so")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so"
         OLD_RPATH "::::::::::::::::::::::::::::::::::::::::::"
         NEW_RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libdestimator.so")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Runtime")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator"
         RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib")
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE FILES "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/test_destimator")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator"
         OLD_RPATH "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src:"
         NEW_RPATH "$ORIGIN/.:$ORIGIN/../lib:$ORIGIN/../../lib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/test_destimator")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  include("/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/FixBundle.cmake")
endif()

