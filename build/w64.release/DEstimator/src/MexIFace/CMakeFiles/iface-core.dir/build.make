# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.2

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator

# Include any dependencies generated for this target.
include src/MexIFace/CMakeFiles/iface-core.dir/depend.make

# Include the progress variables for this target.
include src/MexIFace/CMakeFiles/iface-core.dir/progress.make

# Include the compile flags for this target's objects.
include src/MexIFace/CMakeFiles/iface-core.dir/flags.make

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/flags.make
src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/includes_CXX.rsp
src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/Mex_Iface.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/Mex_Iface.cpp

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/iface-core.dir/Mex_Iface.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/Mex_Iface.cpp > CMakeFiles/iface-core.dir/Mex_Iface.cpp.i

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/iface-core.dir/Mex_Iface.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/Mex_Iface.cpp -o CMakeFiles/iface-core.dir/Mex_Iface.cpp.s

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.requires:
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.requires

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.provides: src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.requires
	$(MAKE) -f src/MexIFace/CMakeFiles/iface-core.dir/build.make src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.provides.build
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.provides

src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.provides.build: src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/flags.make
src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/includes_CXX.rsp
src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/MexUtils.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/iface-core.dir/MexUtils.cpp.obj -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/MexUtils.cpp

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/iface-core.dir/MexUtils.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/MexUtils.cpp > CMakeFiles/iface-core.dir/MexUtils.cpp.i

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/iface-core.dir/MexUtils.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/MexUtils.cpp -o CMakeFiles/iface-core.dir/MexUtils.cpp.s

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.requires:
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.requires

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.provides: src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.requires
	$(MAKE) -f src/MexIFace/CMakeFiles/iface-core.dir/build.make src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.provides.build
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.provides

src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.provides.build: src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/flags.make
src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj: src/MexIFace/CMakeFiles/iface-core.dir/includes_CXX.rsp
src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/explore.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/iface-core.dir/explore.cpp.obj -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/explore.cpp

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/iface-core.dir/explore.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/explore.cpp > CMakeFiles/iface-core.dir/explore.cpp.i

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/iface-core.dir/explore.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && /home/mjo/mxe/usr/bin/x86_64-w64-mingw32.shared-g++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src/explore.cpp -o CMakeFiles/iface-core.dir/explore.cpp.s

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.requires:
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.requires

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.provides: src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.requires
	$(MAKE) -f src/MexIFace/CMakeFiles/iface-core.dir/build.make src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.provides.build
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.provides

src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.provides.build: src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj

iface-core: src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj
iface-core: src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj
iface-core: src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj
iface-core: src/MexIFace/CMakeFiles/iface-core.dir/build.make
.PHONY : iface-core

# Rule to build all files generated by this target.
src/MexIFace/CMakeFiles/iface-core.dir/build: iface-core
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/build

src/MexIFace/CMakeFiles/iface-core.dir/requires: src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.obj.requires
src/MexIFace/CMakeFiles/iface-core.dir/requires: src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.obj.requires
src/MexIFace/CMakeFiles/iface-core.dir/requires: src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.obj.requires
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/requires

src/MexIFace/CMakeFiles/iface-core.dir/clean:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace && $(CMAKE_COMMAND) -P CMakeFiles/iface-core.dir/cmake_clean.cmake
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/clean

src/MexIFace/CMakeFiles/iface-core.dir/depend:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/MexIface/src /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/w64.release/DEstimator/src/MexIFace/CMakeFiles/iface-core.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/MexIFace/CMakeFiles/iface-core.dir/depend

