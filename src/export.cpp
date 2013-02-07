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

  vtkSmartPointer<vtkIntArray> spheresType = vtkSmartPointer<vtkIntArray>::New();
  spheresType->SetNumberOfComponents(1);
  spheresType->SetName("type");
  
  vtkSmartPointer<vtkDoubleArray> spheresVelL = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelL->SetNumberOfComponents(3);
  spheresVelL->SetName("velocity_lin");

  vtkSmartPointer<vtkDoubleArray> spheresVelA = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelA->SetNumberOfComponents(3);
  spheresVelA->SetName("velocity_ang");
  
  vtkSmartPointer<vtkDoubleArray> vectorDr = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDr->SetNumberOfComponents(3);
  vectorDr->SetName("vector_dr");
  
  vtkSmartPointer<vtkDoubleArray> vectorDz = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDz->SetNumberOfComponents(3);
  vectorDz->SetName("vector_dz");
  
  vtkSmartPointer<vtkDoubleArray> vectorDf = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDf->SetNumberOfComponents(3);
  vectorDf->SetName("vector_df");
  
  
  
  for (int z = 0; z<_particleAll->arraySize(); z++) {
    if (_particleAll->particleReal(z)) {
      boost::shared_ptr<particle> partTemp = _particleAll->getP(z);
      vtkIdType pid[1];
      pid[0] = spheresPos->InsertNextPoint(partTemp->c()[0], partTemp->c()[1], partTemp->c()[2]);
      radii->InsertNextValue(partTemp->rad());
      dist->InsertNextValue(partTemp->dist());
      spheresId->InsertNextValue(partTemp->id());
      spheresType->InsertNextValue(partTemp->type());
      
      double vv[3] = {partTemp->v()[0], partTemp->v()[1], partTemp->v()[2]};
      spheresVelL->InsertNextTupleValue(vv);
      
      double aa[3] = {partTemp->o()[0], partTemp->o()[1], partTemp->o()[2]};
      spheresVelA->InsertNextTupleValue(aa);
      
      double dr[3] = {partTemp->dr()[0], partTemp->dr()[1], partTemp->dr()[2]};
      vectorDr->InsertNextTupleValue(dr);
      
      double dz[3] = {partTemp->dz()[0], partTemp->dz()[1], partTemp->dz()[2]};
      vectorDz->InsertNextTupleValue(dz);
      
      double df[3] = {partTemp->df()[0], partTemp->df()[1], partTemp->df()[2]};
      vectorDf->InsertNextTupleValue(df);
      
      
      
      spheresCells->InsertNextCell(1,pid);
    }
  }
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  spheresUg->SetPoints(spheresPos);
  spheresUg->SetCells(VTK_VERTEX, spheresCells);
  spheresUg->GetPointData()->AddArray(radii);
  spheresUg->GetPointData()->AddArray(dist);
  spheresUg->GetPointData()->AddArray(spheresId);
  spheresUg->GetPointData()->AddArray(spheresType);
  spheresUg->GetPointData()->AddArray(spheresVelL);
  spheresUg->GetPointData()->AddArray(spheresVelA);
  spheresUg->GetPointData()->AddArray(vectorDr);
  spheresUg->GetPointData()->AddArray(vectorDz);
  spheresUg->GetPointData()->AddArray(vectorDf);
  
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetInput(spheresUg);
  writer->SetFileName("a.vtu");
  writer->Write();
  
};
