# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


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
CMAKE_SOURCE_DIR = /home/nate/Documents/Repos/LieAlgebra

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/nate/Documents/Repos/LieAlgebra

# Include any dependencies generated for this target.
include libs/CMakeFiles/libs.dir/depend.make

# Include the progress variables for this target.
include libs/CMakeFiles/libs.dir/progress.make

# Include the compile flags for this target's objects.
include libs/CMakeFiles/libs.dir/flags.make

libs/CMakeFiles/libs.dir/RootStructure.cpp.o: libs/CMakeFiles/libs.dir/flags.make
libs/CMakeFiles/libs.dir/RootStructure.cpp.o: libs/RootStructure.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/nate/Documents/Repos/LieAlgebra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object libs/CMakeFiles/libs.dir/RootStructure.cpp.o"
	cd /home/nate/Documents/Repos/LieAlgebra/libs && /usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/libs.dir/RootStructure.cpp.o -c /home/nate/Documents/Repos/LieAlgebra/libs/RootStructure.cpp

libs/CMakeFiles/libs.dir/RootStructure.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/libs.dir/RootStructure.cpp.i"
	cd /home/nate/Documents/Repos/LieAlgebra/libs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/nate/Documents/Repos/LieAlgebra/libs/RootStructure.cpp > CMakeFiles/libs.dir/RootStructure.cpp.i

libs/CMakeFiles/libs.dir/RootStructure.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/libs.dir/RootStructure.cpp.s"
	cd /home/nate/Documents/Repos/LieAlgebra/libs && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/nate/Documents/Repos/LieAlgebra/libs/RootStructure.cpp -o CMakeFiles/libs.dir/RootStructure.cpp.s

libs/CMakeFiles/libs.dir/RootStructure.cpp.o.requires:

.PHONY : libs/CMakeFiles/libs.dir/RootStructure.cpp.o.requires

libs/CMakeFiles/libs.dir/RootStructure.cpp.o.provides: libs/CMakeFiles/libs.dir/RootStructure.cpp.o.requires
	$(MAKE) -f libs/CMakeFiles/libs.dir/build.make libs/CMakeFiles/libs.dir/RootStructure.cpp.o.provides.build
.PHONY : libs/CMakeFiles/libs.dir/RootStructure.cpp.o.provides

libs/CMakeFiles/libs.dir/RootStructure.cpp.o.provides.build: libs/CMakeFiles/libs.dir/RootStructure.cpp.o


# Object files for target libs
libs_OBJECTS = \
"CMakeFiles/libs.dir/RootStructure.cpp.o"

# External object files for target libs
libs_EXTERNAL_OBJECTS =

libs/liblibs.a: libs/CMakeFiles/libs.dir/RootStructure.cpp.o
libs/liblibs.a: libs/CMakeFiles/libs.dir/build.make
libs/liblibs.a: libs/CMakeFiles/libs.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/nate/Documents/Repos/LieAlgebra/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library liblibs.a"
	cd /home/nate/Documents/Repos/LieAlgebra/libs && $(CMAKE_COMMAND) -P CMakeFiles/libs.dir/cmake_clean_target.cmake
	cd /home/nate/Documents/Repos/LieAlgebra/libs && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/libs.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
libs/CMakeFiles/libs.dir/build: libs/liblibs.a

.PHONY : libs/CMakeFiles/libs.dir/build

libs/CMakeFiles/libs.dir/requires: libs/CMakeFiles/libs.dir/RootStructure.cpp.o.requires

.PHONY : libs/CMakeFiles/libs.dir/requires

libs/CMakeFiles/libs.dir/clean:
	cd /home/nate/Documents/Repos/LieAlgebra/libs && $(CMAKE_COMMAND) -P CMakeFiles/libs.dir/cmake_clean.cmake
.PHONY : libs/CMakeFiles/libs.dir/clean

libs/CMakeFiles/libs.dir/depend:
	cd /home/nate/Documents/Repos/LieAlgebra && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/nate/Documents/Repos/LieAlgebra /home/nate/Documents/Repos/LieAlgebra/libs /home/nate/Documents/Repos/LieAlgebra /home/nate/Documents/Repos/LieAlgebra/libs /home/nate/Documents/Repos/LieAlgebra/libs/CMakeFiles/libs.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : libs/CMakeFiles/libs.dir/depend

