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
CMAKE_BINARY_DIR = /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator

# Include any dependencies generated for this target.
include src/CMakeFiles/DEstimator_Iface.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/DEstimator_Iface.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/DEstimator_Iface.dir/flags.make

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o: src/CMakeFiles/DEstimator_Iface.dir/flags.make
src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DEstimator_Iface.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DEstimator_Iface.cpp

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DEstimator_Iface.cpp > CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.i

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DEstimator_Iface.cpp -o CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.s

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.requires:
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.requires

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.provides: src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/DEstimator_Iface.dir/build.make src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.provides.build
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.provides

src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.provides.build: src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o

# Object files for target DEstimator_Iface
DEstimator_Iface_OBJECTS = \
"CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o"

# External object files for target DEstimator_Iface
DEstimator_Iface_EXTERNAL_OBJECTS = \
"/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.o" \
"/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.o" \
"/home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.o"

src/DEstimator_Iface.mexa64: src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o
src/DEstimator_Iface.mexa64: src/MexIFace/CMakeFiles/iface-core.dir/Mex_Iface.cpp.o
src/DEstimator_Iface.mexa64: src/MexIFace/CMakeFiles/iface-core.dir/MexUtils.cpp.o
src/DEstimator_Iface.mexa64: src/MexIFace/CMakeFiles/iface-core.dir/explore.cpp.o
src/DEstimator_Iface.mexa64: src/CMakeFiles/DEstimator_Iface.dir/build.make
src/DEstimator_Iface.mexa64: /usr/lib64/librefblas.so
src/DEstimator_Iface.mexa64: /opt/MATLAB/R2014b/bin/glnxa64/libmex.so
src/DEstimator_Iface.mexa64: /opt/MATLAB/R2014b/bin/glnxa64/libmx.so
src/DEstimator_Iface.mexa64: /opt/MATLAB/R2014b/bin/glnxa64/libeng.so
src/DEstimator_Iface.mexa64: /opt/MATLAB/R2014b/bin/glnxa64/libmat.so
src/DEstimator_Iface.mexa64: /usr/lib64/libpthread.so
src/DEstimator_Iface.mexa64: src/libdestimator.so
src/DEstimator_Iface.mexa64: /usr/lib64/librefblas.so
src/DEstimator_Iface.mexa64: src/CMakeFiles/DEstimator_Iface.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library DEstimator_Iface.mexa64"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/DEstimator_Iface.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/DEstimator_Iface.dir/build: src/DEstimator_Iface.mexa64
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/build

src/CMakeFiles/DEstimator_Iface.dir/requires: src/CMakeFiles/DEstimator_Iface.dir/DEstimator_Iface.cpp.o.requires
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/requires

src/CMakeFiles/DEstimator_Iface.dir/clean:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && $(CMAKE_COMMAND) -P CMakeFiles/DEstimator_Iface.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/clean

src/CMakeFiles/DEstimator_Iface.dir/depend:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/CMakeFiles/DEstimator_Iface.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/DEstimator_Iface.dir/depend

