# Install script for directory: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "../../../mex/mex.w64")
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

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Runtime")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE SHARED_LIBRARY FILES "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/DEstimator_Iface.mexw64")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./DEstimator_Iface.mexw64" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./DEstimator_Iface.mexw64")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./DEstimator_Iface.mexw64")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Runtime")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE SHARED_LIBRARY FILES "/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/libdestimator.dll")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./libdestimator.dll" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./libdestimator.dll")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-strip" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/./libdestimator.dll")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  include("/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/FixBundle.cmake")
endif()

