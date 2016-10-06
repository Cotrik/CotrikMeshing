# CMake entry point
cmake_minimum_required (VERSION 2.6)
project (Hexmeshing)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenGL REQUIRED)
find_package(Boost QUIET COMPONENTS program_options)
if(Boost_PROGRAM_OPTIONS_FOUND)
    if( CGAL_AUTO_LINK_ENABLED )
        message( STATUS "Boost.ProgramOptions library: found" )
    else()
        message( STATUS "Boost.ProgramOptions library: ${Boost_PROGRAM_OPTIONS_LIBRARY}" )
    endif()
    add_definitions( "-DCGAL_USE_BOOST_PROGRAM_OPTIONS" )
    list(APPEND CGAL_3RD_PARTY_LIBRARIES ${Boost_LIBRARIES})
endif()

#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

if( CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR )
    message( FATAL_ERROR "Please select another Build Directory ! (and give it a clever name, like bin_Visual2012_64bits/)" )
endif()
if( CMAKE_SOURCE_DIR MATCHES " " )
	message( "Your Source Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()
if( CMAKE_BINARY_DIR MATCHES " " )
	message( "Your Build Directory contains spaces. If you experience problems when compiling, this can be the cause." )
endif()

find_package(CGAL QUIET COMPONENTS Core )

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )

  include_directories (BEFORE "../../include")
endif()
# Compile external dependencies
# add_subdirectory (external)
# add_subdirectory( GeodesicsOnMeshes/FastMarching )
# add_subdirectory( GeodesicsOnMeshes/Source )
# add_subdirectory( GeodesicsOnMeshes/Examples )
# add_subdirectory( GeodesicsOnMeshes/Testing )

# On Visual 2005 and above, this module can set the debug working directory
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/rpavlik-cmake-modules-1c73e35")
#include(CreateLaunchers)
#include(MSVCMultipleProcessCompile) # /MP

if(INCLUDE_DISTRIB)
	add_subdirectory(distrib)
endif(INCLUDE_DISTRIB)

include_directories(
#	external/glfw-3.0.3/include/GLFW/
#	external/glew-1.9.0/include/
    /usr/include/
    /usr/include/eigen3
    /usr/local/include/
#    /usr/local/include/vtk-6.3/
    src/
	.
)

set(ALL_LIBS
	${OPENGL_LIBRARY}
#	GLFW_303
#	GLEW_190
#	ANTTWEAKBAR_116_OGLCORE_GLFW
	GL
    glut
    boost_program_options
    pthread
#     ${VTK_LIBRARIES}
#     vtkMeshGeodesics
)

add_definitions(
#	-DTW_STATIC
#	-DTW_NO_LIB_PRAGMA
#	-DTW_NO_DIRECT3D
#	-DGLEW_STATIC
#	-D_CRT_SECURE_NO_WARNINGS
)

set(SRC
    src/include.cpp
    src/include.h
#	src/FocusContextMagnifier.cpp
#    src/FocusContextMagnifier.h
#    src/SimpleGeometryViewer.h
#    src/SimpleGeometryViewer.cpp
    src/GeometricStruct.cpp
    src/GeometricStruct.h
    src/GeoUtil.cpp
	src/GeoUtil.h
#	src/Iterator.h
#	src/Log.cpp
#	src/Log.h
	src/Mesh.cpp
#	src/MeshDisplayer.cpp
#	src/MeshDisplayer.h
	src/MeshFileReader.cpp
	src/MeshFileReader.h
	src/MeshFileWriter.cpp
	src/MeshFileWriter.h
	src/Mesh.h
	src/PolycubeMesh.cpp
	src/PolycubeMesh.h
#	src/math3d.cpp
#	src/math3d.h
#	src/Selection.cpp
#	src/Selection.h
	src/Patch.cpp
	src/Patch.h
	src/Parametrizer.cpp
	src/Parametrizer.h
)

#add_executable(tetdensity ${SRC} app/tetdensity.cpp)
#target_link_libraries(tetdensity	${ALL_LIBS})
#add_executable(hexdensity ${SRC} app/hexdensity.cpp)
#target_link_libraries(hexdensity	${ALL_LIBS})
#add_executable(magnify ${SRC}	app/magnify.cpp)
#target_link_libraries(magnify	${ALL_LIBS})
#add_executable(magnify_Cp ${SRC}	app/magnify_Cp.cpp)
#target_link_libraries(magnify_Cp	${ALL_LIBS})
#add_executable(meshdilate ${SRC}	app/meshdilate.cpp)
#target_link_libraries(meshdilate	${ALL_LIBS})
#add_executable(showedge ${SRC} app/showedge.cpp)
#target_link_libraries(showedge	${ALL_LIBS})
#add_executable(test_llt_qr_lu ${SRC} app/testLLT_QP_LU.cpp)
#target_link_libraries(test_llt_qr_lu ${ALL_LIBS})
#add_executable(polycube_cleanup ${SRC} app/polycube_cleanup.cpp)
#target_link_libraries(polycube_cleanup ${ALL_LIBS})
#add_executable(hexgen_plus ${SRC} app/hexgen_plus.cpp)
#target_link_libraries(hexgen_plus ${ALL_LIBS})
#add_executable(hexgen_fast ${SRC} app/hexgen_fast.cpp)
#target_link_libraries(hexgen_fast ${ALL_LIBS})
#add_executable(polycube_hexmeshing ${SRC} app/polycube_hexmeshing.cpp)
#target_link_libraries(polycube_hexmeshing ${ALL_LIBS})
#add_executable(display_mesh ${SRC} app/display_mesh.cpp)
#target_link_libraries(display_mesh ${ALL_LIBS})
#add_executable(display_highlight_patch ${SRC} app/display_highlight_patch.cpp)
#target_link_libraries(display_highlight_patch ${ALL_LIBS})
#add_executable(detect_overlapping_area ${SRC} app/detect_overlapping_area.cpp)
#target_link_libraries(detect_overlapping_area ${ALL_LIBS})
#add_executable(ajust_polycube ${SRC} app/ajust_polycube.cpp)
#target_link_libraries(ajust_polycube ${ALL_LIBS})

#add_executable(calc_scaled_jacobian	${SRC}	app/calc_scaled_jacobian.cpp)
#target_link_libraries(calc_scaled_jacobian	${ALL_LIBS})
#add_executable(hexmapback ${SRC} app/hexmapback.cpp)
#target_link_libraries(hexmapback ${ALL_LIBS})
#add_executable(hexgen ${SRC} app/hexgen.cpp)
#target_link_libraries(hexgen ${ALL_LIBS})
#add_executable(hexgen_align ${SRC} app/hexgen_align.cpp)
#target_link_libraries(hexgen_align ${ALL_LIBS})
#add_executable(hexgen_auto ${SRC} app/hexgen_auto.cpp)
#target_link_libraries(hexgen_auto ${ALL_LIBS})
#add_executable(hexsmooth ${SRC} app/hexsmooth.cpp)
#target_link_libraries(hexsmooth ${ALL_LIBS})
#add_executable(vol2vtk ${SRC} app/vol2vtk.cpp)
#target_link_libraries(vol2vtk ${ALL_LIBS})
#add_executable(align_polycube_surface ${SRC} app/align_polycube_surface.cpp)
#target_link_libraries(align_polycube_surface ${ALL_LIBS})
#add_executable(align_polycube_surface_patch ${SRC} app/align_polycube_surface_patch.cpp)
#target_link_libraries(align_polycube_surface_patch ${ALL_LIBS})
#add_executable(align_polycube_surface_patch2 ${SRC} app/align_polycube_surface_patch2.cpp)
#target_link_libraries(align_polycube_surface_patch2 ${ALL_LIBS})
add_executable(convert_file ${SRC} app/convert_file.cpp)
target_link_libraries(convert_file ${ALL_LIBS})
#add_executable(verify_structure ${SRC} app/verify_structure.cpp)
#target_link_libraries(verify_structure ${ALL_LIBS})
#add_executable(reorient ${SRC} app/reorient.cpp)
#target_link_libraries(reorient ${ALL_LIBS})
#add_executable(cone_opt ${SRC} app/cone_optimize.cpp)
#target_link_libraries(cone_opt ${ALL_LIBS})
#add_executable(gen_rot_meshes ${SRC} app/gen_rot_meshes.cpp)
#target_link_libraries(gen_rot_meshes ${ALL_LIBS})
#add_executable(hex_refine ${SRC} app/hex_refine.cpp)
#target_link_libraries(hex_refine ${ALL_LIBS})
#add_executable(boundary ${SRC} app/boundary.cpp)
#target_link_libraries(boundary ${ALL_LIBS})
#add_executable(normalize ${SRC} app/normalize.cpp)
#target_link_libraries(normalize ${ALL_LIBS})
#add_executable(extract_surface ${SRC} app/extract_surface.cpp)
#target_link_libraries(extract_surface ${ALL_LIBS})
#add_executable(extract_singularities ${SRC} app/extract_singularities.cpp)
#target_link_libraries(extract_singularities ${ALL_LIBS})
#add_executable(hexgen_s ${SRC} app/hexgen_s.cpp)
#target_link_libraries(hexgen_s ${ALL_LIBS})
#add_executable(polycubehexgen ${SRC} app/polycubehexgen.cpp)
#target_link_libraries(polycubehexgen ${ALL_LIBS})
#add_executable(polycubehexgen_round ${SRC} app/polycubehexgen_round.cpp)
#target_link_libraries(polycubehexgen_round ${ALL_LIBS})
#add_executable(polycubehexgen_round2 ${SRC} app/polycubehexgen_round2.cpp)
#target_link_libraries(polycubehexgen_round2 ${ALL_LIBS})
#add_executable(polycubehexgen_aniso ${SRC} app/polycubehexgen_aniso.cpp)
#target_link_libraries(polycubehexgen_aniso ${ALL_LIBS})
#add_executable(polycubehexgen_round_aniso ${SRC} app/polycubehexgen_round_aniso.cpp)
#target_link_libraries(polycubehexgen_round_aniso ${ALL_LIBS})

#add_executable(polycubehexgen_round_arap_aniso ${SRC} app/polycubehexgen_round_arap_aniso.cpp)
#target_link_libraries(polycubehexgen_round_arap_aniso ${ALL_LIBS})

#add_executable(hexgen_surface ${SRC} app/hexgen_surface.cpp)
#target_link_libraries(hexgen_surface ${ALL_LIBS})
#add_executable(magnify_segments ${SRC} app/magnify_segments.cpp)
#target_link_libraries(magnify_segments ${ALL_LIBS})
#add_executable(magnify_multi_segments ${SRC} app/magnify_multi_segments.cpp)
#target_link_libraries(magnify_multi_segments ${ALL_LIBS})
#add_executable(magnify_multi_segments_auto ${SRC} app/magnify_multi_segments_auto.cpp)
#target_link_libraries(magnify_multi_segments_auto ${ALL_LIBS})
#add_executable(magnify_multi_segments_scales ${SRC} app/magnify_multi_segments_scales.cpp)
#target_link_libraries(magnify_multi_segments_scales ${ALL_LIBS})

#add_executable(magnify_multi_segments_scales_aniso ${SRC} app/magnify_multi_segments_scales_aniso.cpp)
#target_link_libraries(magnify_multi_segments_scales_aniso ${ALL_LIBS})

#add_executable(magnify_multi_segments_scales_volume_aniso ${SRC} app/magnify_multi_segments_scales_volume_aniso.cpp)
#target_link_libraries(magnify_multi_segments_scales_volume_aniso ${ALL_LIBS})

#add_executable(magnify_multi_segments_scales_volume_aniso2 ${SRC} app/magnify_multi_segments_scales_volume_aniso2.cpp)
#target_link_libraries(magnify_multi_segments_scales_volume_aniso2 ${ALL_LIBS})

#add_executable(magnify_multi_segments_scales_volume_aniso3 ${SRC} app/magnify_multi_segments_scales_volume_aniso3.cpp)
#target_link_libraries(magnify_multi_segments_scales_volume_aniso3 ${ALL_LIBS})

#add_executable(magnify_segments_test ${SRC} app/magnify_segments_test.cpp)
#target_link_libraries(magnify_segments_test ${ALL_LIBS})
#add_executable(mesh_map ${SRC} app/mesh_map.cpp)
#target_link_libraries(mesh_map ${ALL_LIBS})
#add_executable(hex_map ${SRC} app/hex_map.cpp)
#target_link_libraries(hex_map ${ALL_LIBS})
#add_executable(check_overlap ${SRC} app/check_overlap.cpp)
#target_link_libraries(check_overlap ${ALL_LIBS})
# add_executable(check_overlap_clean ${SRC} app/check_overlap_clean.cpp)
# target_link_libraries(check_overlap_clean ${ALL_LIBS})
#add_executable(hexgen_clean ${SRC} app/hexgen_clean.cpp)
#target_link_libraries(hexgen_clean ${ALL_LIBS})
#add_executable(arap_tet ${SRC} app/arap_tet.cpp)
#target_link_libraries(arap_tet ${ALL_LIBS})
#add_executable(magnify_segments_volume ${SRC} app/magnify_segments_volume.cpp)
#target_link_libraries(magnify_segments_volume ${ALL_LIBS})
#add_executable(mean_value_coordinates_parametrization ${SRC} app/mean_value_coordinates_parametrization.cpp)
#target_link_libraries(mean_value_coordinates_parametrization ${ALL_LIBS})

# add_executable(polycube_hex_extract ${SRC} app/polycube_hex_extract.cpp)
# target_link_libraries(polycube_hex_extract ${ALL_LIBS})