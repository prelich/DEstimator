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
include src/CMakeFiles/destimator.dir/depend.make

# Include the progress variables for this target.
include src/CMakeFiles/destimator.dir/progress.make

# Include the compile flags for this target's objects.
include src/CMakeFiles/destimator.dir/flags.make

src/CMakeFiles/destimator.dir/destimator.cpp.o: src/CMakeFiles/destimator.dir/flags.make
src/CMakeFiles/destimator.dir/destimator.cpp.o: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/destimator.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/destimator.dir/destimator.cpp.o"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/destimator.dir/destimator.cpp.o -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/destimator.cpp

src/CMakeFiles/destimator.dir/destimator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/destimator.dir/destimator.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/destimator.cpp > CMakeFiles/destimator.dir/destimator.cpp.i

src/CMakeFiles/destimator.dir/destimator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/destimator.dir/destimator.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/destimator.cpp -o CMakeFiles/destimator.dir/destimator.cpp.s

src/CMakeFiles/destimator.dir/destimator.cpp.o.requires:
.PHONY : src/CMakeFiles/destimator.dir/destimator.cpp.o.requires

src/CMakeFiles/destimator.dir/destimator.cpp.o.provides: src/CMakeFiles/destimator.dir/destimator.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/destimator.dir/build.make src/CMakeFiles/destimator.dir/destimator.cpp.o.provides.build
.PHONY : src/CMakeFiles/destimator.dir/destimator.cpp.o.provides

src/CMakeFiles/destimator.dir/destimator.cpp.o.provides.build: src/CMakeFiles/destimator.dir/destimator.cpp.o

src/CMakeFiles/destimator.dir/DCompLib.cpp.o: src/CMakeFiles/destimator.dir/flags.make
src/CMakeFiles/destimator.dir/DCompLib.cpp.o: /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DCompLib.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/CMakeFiles/destimator.dir/DCompLib.cpp.o"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/destimator.dir/DCompLib.cpp.o -c /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DCompLib.cpp

src/CMakeFiles/destimator.dir/DCompLib.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/destimator.dir/DCompLib.cpp.i"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DCompLib.cpp > CMakeFiles/destimator.dir/DCompLib.cpp.i

src/CMakeFiles/destimator.dir/DCompLib.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/destimator.dir/DCompLib.cpp.s"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src/DCompLib.cpp -o CMakeFiles/destimator.dir/DCompLib.cpp.s

src/CMakeFiles/destimator.dir/DCompLib.cpp.o.requires:
.PHONY : src/CMakeFiles/destimator.dir/DCompLib.cpp.o.requires

src/CMakeFiles/destimator.dir/DCompLib.cpp.o.provides: src/CMakeFiles/destimator.dir/DCompLib.cpp.o.requires
	$(MAKE) -f src/CMakeFiles/destimator.dir/build.make src/CMakeFiles/destimator.dir/DCompLib.cpp.o.provides.build
.PHONY : src/CMakeFiles/destimator.dir/DCompLib.cpp.o.provides

src/CMakeFiles/destimator.dir/DCompLib.cpp.o.provides.build: src/CMakeFiles/destimator.dir/DCompLib.cpp.o

# Object files for target destimator
destimator_OBJECTS = \
"CMakeFiles/destimator.dir/destimator.cpp.o" \
"CMakeFiles/destimator.dir/DCompLib.cpp.o"

# External object files for target destimator
destimator_EXTERNAL_OBJECTS =

src/libdestimator.so: src/CMakeFiles/destimator.dir/destimator.cpp.o
src/libdestimator.so: src/CMakeFiles/destimator.dir/DCompLib.cpp.o
src/libdestimator.so: src/CMakeFiles/destimator.dir/build.make
src/libdestimator.so: /usr/lib64/librefblas.so
src/libdestimator.so: src/CMakeFiles/destimator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared library libdestimator.so"
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/destimator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/CMakeFiles/destimator.dir/build: src/libdestimator.so
.PHONY : src/CMakeFiles/destimator.dir/build

src/CMakeFiles/destimator.dir/requires: src/CMakeFiles/destimator.dir/destimator.cpp.o.requires
src/CMakeFiles/destimator.dir/requires: src/CMakeFiles/destimator.dir/DCompLib.cpp.o.requires
.PHONY : src/CMakeFiles/destimator.dir/requires

src/CMakeFiles/destimator.dir/clean:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src && $(CMAKE_COMMAND) -P CMakeFiles/destimator.dir/cmake_clean.cmake
.PHONY : src/CMakeFiles/destimator.dir/clean

src/CMakeFiles/destimator.dir/depend:
	cd /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/tools/DEstimator/src /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src /home/mjo/temp/SMA_Toolbox/tags/DEstimator-1.0/build/linux.release/DEstimator/src/CMakeFiles/destimator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/CMakeFiles/destimator.dir/depend

