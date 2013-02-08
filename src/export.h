#ifndef EXPORTCLASS
#define EXPORTCLASS

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>


#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkTriangle.h>
#include <vtkLine.h>
#include <vtkQuad.h>


#include "config.h"
#include "particle.h"
#include "force.h"
#include "band.h"



using namespace std;
class exportclass {
  private:
    boost::shared_ptr <configopt> _cfg;
    boost::shared_ptr <bandRow> _bandRow;
    ofstream fileVTK;
    const char * _fileName;
    
  public:
    exportclass(boost::shared_ptr<configopt>, boost::shared_ptr <bandRow>);
    void exportVTK();
};

#endif
