#ifndef __INCLUDE_H__
#define __INCLUDE_H__
#include "Iterator.h"
#include "Mesh.h"
#include "MeshFileReader.h"
#include "MeshFileWriter.h"
#include "MeshDisplayer.h"
#include "PolycubeMesh.h"
#include "FocusContextMagnifier.h"
#include "SimpleGeometryViewer.h"
#include <sstream>
#include <iostream>
#include <exception>
#include <GL/glut.h>
#ifdef CGAL_USE_BOOST_PROGRAM_OPTIONS
#include <boost/program_options.hpp>
namespace po = boost::program_options;
#endif

#ifdef _WIN32
#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Sparse"
//#include "Eigen/MPRealSupport"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/SparseCore"
//#include "Eigen/src/MPRealSupport/mpfr/mpfr.h"
#else
#include "Eigen/Core"
#include "Eigen/Eigen"
#include "Eigen/Dense"
#include "Eigen/Sparse"
//#include "eigen3/Eigen/MPRealSupport"
#include "Eigen/Cholesky"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Eigen/SparseCore"
#endif


#include <vtkCellType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCellIterator.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyData.h>

#include "vtkSmartPointer.h"

#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkTimerLog.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPolyDataNormals.h"
#include "vtkRendererCollection.h"
#include "vtkPolyDataCollection.h"
#include "vtkObjectFactory.h"
#include "vtkIdList.h"
#include "vtkPolyDataReader.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkNew.h"
#include "vtkPointData.h"
#include "vtkContourWidget.h"
#include "vtkOrientedGlyphContourRepresentation.h"
#include "vtkPolygonalSurfacePointPlacer.h"
#include "../GeodesicsOnMeshes/Source/vtkPolygonalSurfaceContourLineInterpolator2.h"
extern float CUBESIZE;

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
using namespace std;

class Timer {
private:

    timeval startTime;

public:

    void start(){
        gettimeofday(&startTime, NULL);
    }

    double stop(){
        timeval endTime;
        long seconds, useconds;
        double duration;

        gettimeofday(&endTime, NULL);

        seconds  = endTime.tv_sec  - startTime.tv_sec;
        useconds = endTime.tv_usec - startTime.tv_usec;

        duration = seconds + useconds/1000000.0;

        return duration;
    }

    void printTime(double duration){
        printf("%5.6f seconds\n", duration);
    }
};
#endif // __INCLUDE_H__
