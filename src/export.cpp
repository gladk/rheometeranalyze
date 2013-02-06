#include "export.h"

exportclass::exportclass(boost::shared_ptr<configopt> cfg, boost::shared_ptr <particleRow> particleAll) {
  _cfg = cfg;
  _particleAll = particleAll;
};

void exportclass::exportVTK() {
  vtkSmartPointer<vtkPoints>  spheresPos = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("radii");

  vtkSmartPointer<vtkDoubleArray> dist = vtkSmartPointer<vtkDoubleArray>::New();
  dist->SetNumberOfComponents(1);
  dist->SetName("dist");
  
  vtkSmartPointer<vtkIntArray> spheresId = vtkSmartPointer<vtkIntArray>::New();
  spheresId->SetNumberOfComponents(1);
  spheresId->SetName("id");
  
  for (int z = 0; z<_particleAll->arraySize(); z++) {
    if (_particleAll->particleReal(z)) {
      boost::shared_ptr<particle> partTemp = _particleAll->getP(z);
      vtkIdType pid[1];
      pid[0] = spheresPos->InsertNextPoint(partTemp->c()[0], partTemp->c()[1], partTemp->c()[2]);
      radii->InsertNextValue(partTemp->rad());
      dist->InsertNextValue(partTemp->dist());
      spheresId->InsertNextValue(partTemp->id());
      
      spheresCells->InsertNextCell(1,pid);
    }
  }
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  spheresUg->SetPoints(spheresPos);
  spheresUg->SetCells(VTK_VERTEX, spheresCells);
  spheresUg->GetPointData()->AddArray(radii);
  spheresUg->GetPointData()->AddArray(dist);
  spheresUg->GetPointData()->AddArray(spheresId);
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetInput(spheresUg);
  writer->SetFileName("a.vtk");
  writer->Write();
  
};
