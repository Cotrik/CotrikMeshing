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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/cotrik/Downloads/GeodesicsOnMeshes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/cotrik/Downloads/GeodesicsOnMeshes/build

# Include any dependencies generated for this target.
include Source/CMakeFiles/vtkMeshGeodesics.dir/depend.make

# Include the progress variables for this target.
include Source/CMakeFiles/vtkMeshGeodesics.dir/progress.make

# Include the compile flags for this target's objects.
include Source/CMakeFiles/vtkMeshGeodesics.dir/flags.make

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o: Source/CMakeFiles/vtkMeshGeodesics.dir/flags.make
Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o: ../Source/vtkPolyDataGeodesicDistance.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cotrik/Downloads/GeodesicsOnMeshes/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o -c /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolyDataGeodesicDistance.cxx

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.i"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolyDataGeodesicDistance.cxx > CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.i

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.s"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolyDataGeodesicDistance.cxx -o CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.s

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.requires:
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.requires

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.provides: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.requires
	$(MAKE) -f Source/CMakeFiles/vtkMeshGeodesics.dir/build.make Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.provides.build
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.provides

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.provides.build: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o: Source/CMakeFiles/vtkMeshGeodesics.dir/flags.make
Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o: ../Source/vtkFastMarchingGeodesicDistance.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cotrik/Downloads/GeodesicsOnMeshes/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o -c /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicDistance.cxx

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.i"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicDistance.cxx > CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.i

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.s"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicDistance.cxx -o CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.s

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.requires:
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.requires

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.provides: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.requires
	$(MAKE) -f Source/CMakeFiles/vtkMeshGeodesics.dir/build.make Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.provides.build
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.provides

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.provides.build: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o: Source/CMakeFiles/vtkMeshGeodesics.dir/flags.make
Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o: ../Source/vtkFastMarchingGeodesicPath.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cotrik/Downloads/GeodesicsOnMeshes/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o -c /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicPath.cxx

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.i"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicPath.cxx > CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.i

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.s"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkFastMarchingGeodesicPath.cxx -o CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.s

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.requires:
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.requires

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.provides: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.requires
	$(MAKE) -f Source/CMakeFiles/vtkMeshGeodesics.dir/build.make Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.provides.build
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.provides

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.provides.build: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o: Source/CMakeFiles/vtkMeshGeodesics.dir/flags.make
Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o: ../Source/vtkPolygonalSurfaceContourLineInterpolator2.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /home/cotrik/Downloads/GeodesicsOnMeshes/build/CMakeFiles $(CMAKE_PROGRESS_4)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o -c /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolygonalSurfaceContourLineInterpolator2.cxx

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.i"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolygonalSurfaceContourLineInterpolator2.cxx > CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.i

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.s"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && /usr/bin/c++  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/cotrik/Downloads/GeodesicsOnMeshes/Source/vtkPolygonalSurfaceContourLineInterpolator2.cxx -o CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.s

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.requires:
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.requires

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.provides: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.requires
	$(MAKE) -f Source/CMakeFiles/vtkMeshGeodesics.dir/build.make Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.provides.build
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.provides

Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.provides.build: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o

# Object files for target vtkMeshGeodesics
vtkMeshGeodesics_OBJECTS = \
"CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o" \
"CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o" \
"CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o" \
"CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o"

# External object files for target vtkMeshGeodesics
vtkMeshGeodesics_EXTERNAL_OBJECTS =

bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o
bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o
bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o
bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o
bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/build.make
bin/libvtkMeshGeodesics.a: Source/CMakeFiles/vtkMeshGeodesics.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX static library ../bin/libvtkMeshGeodesics.a"
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && $(CMAKE_COMMAND) -P CMakeFiles/vtkMeshGeodesics.dir/cmake_clean_target.cmake
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vtkMeshGeodesics.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
Source/CMakeFiles/vtkMeshGeodesics.dir/build: bin/libvtkMeshGeodesics.a
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/build

Source/CMakeFiles/vtkMeshGeodesics.dir/requires: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolyDataGeodesicDistance.cxx.o.requires
Source/CMakeFiles/vtkMeshGeodesics.dir/requires: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicDistance.cxx.o.requires
Source/CMakeFiles/vtkMeshGeodesics.dir/requires: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkFastMarchingGeodesicPath.cxx.o.requires
Source/CMakeFiles/vtkMeshGeodesics.dir/requires: Source/CMakeFiles/vtkMeshGeodesics.dir/vtkPolygonalSurfaceContourLineInterpolator2.cxx.o.requires
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/requires

Source/CMakeFiles/vtkMeshGeodesics.dir/clean:
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source && $(CMAKE_COMMAND) -P CMakeFiles/vtkMeshGeodesics.dir/cmake_clean.cmake
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/clean

Source/CMakeFiles/vtkMeshGeodesics.dir/depend:
	cd /home/cotrik/Downloads/GeodesicsOnMeshes/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/cotrik/Downloads/GeodesicsOnMeshes /home/cotrik/Downloads/GeodesicsOnMeshes/Source /home/cotrik/Downloads/GeodesicsOnMeshes/build /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source /home/cotrik/Downloads/GeodesicsOnMeshes/build/Source/CMakeFiles/vtkMeshGeodesics.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : Source/CMakeFiles/vtkMeshGeodesics.dir/depend

