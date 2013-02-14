#ifndef EXPORTCLASS
#define EXPORTCLASS

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <memory>


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
    std::shared_ptr <configopt> _cfg;
    std::shared_ptr <bandRow> _bandRow;
    ofstream fileVTK;
    const char * _fileName;
    
  public:
    exportclass(std::shared_ptr<configopt>, std::shared_ptr <bandRow>);
    void VTK();
    void gnuplotSchearRate();
};

#endif
