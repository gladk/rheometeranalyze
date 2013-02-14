#include "export.h"

exportclass::exportclass(std::shared_ptr<configopt> cfg, std::shared_ptr <bandRow> bandAll) {
  _cfg = cfg;
  _bandRow = bandAll;
};

void exportclass::VTK() {
  
  //Export Particles
  vtkSmartPointer<vtkPoints>  spheresPos = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("radii");

  vtkSmartPointer<vtkDoubleArray> dist = vtkSmartPointer<vtkDoubleArray>::New();
  dist->SetNumberOfComponents(1);
  dist->SetName("dist");

  vtkSmartPointer<vtkDoubleArray> height = vtkSmartPointer<vtkDoubleArray>::New();
  height->SetNumberOfComponents(1);
  height->SetName("height");
  
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
  
  vtkSmartPointer<vtkIntArray> bandR = vtkSmartPointer<vtkIntArray>::New();
  bandR->SetNumberOfComponents(1);
  bandR->SetName("bandR");
  
  vtkSmartPointer<vtkIntArray> bandZ = vtkSmartPointer<vtkIntArray>::New();
  bandZ->SetNumberOfComponents(1);
  bandZ->SetName("bandZ");
  
  vtkSmartPointer<vtkIntArray> bandN = vtkSmartPointer<vtkIntArray>::New();
  bandN->SetNumberOfComponents(1);
  bandN->SetName("bandN");
  
  vtkSmartPointer<vtkDoubleArray> bandTau = vtkSmartPointer<vtkDoubleArray>::New();
  bandTau->SetNumberOfComponents(1);
  bandTau->SetName("bandGlobTau");
  
  vtkSmartPointer<vtkDoubleArray> bandPress = vtkSmartPointer<vtkDoubleArray>::New();
  bandPress->SetNumberOfComponents(1);
  bandPress->SetName("bandGlobPress");
  
  vtkSmartPointer<vtkDoubleArray> bandLocalPress = vtkSmartPointer<vtkDoubleArray>::New();
  bandLocalPress->SetNumberOfComponents(1);
  bandLocalPress->SetName("bandLocPress");
  
  vtkSmartPointer<vtkDoubleArray> bandLocSigmDev = vtkSmartPointer<vtkDoubleArray>::New();
  bandLocSigmDev->SetNumberOfComponents(1);
  bandLocSigmDev->SetName("bandLocSigmDev");
  
  vtkSmartPointer<vtkDoubleArray> bandLocMu = vtkSmartPointer<vtkDoubleArray>::New();
  bandLocMu->SetNumberOfComponents(1);
  bandLocMu->SetName("bandLocMu");
  
  vtkSmartPointer<vtkDoubleArray> bandOmega = vtkSmartPointer<vtkDoubleArray>::New();
  bandOmega->SetNumberOfComponents(1);
  bandOmega->SetName("bandOmega");
  
  vtkSmartPointer<vtkIntArray> bandPartNum = vtkSmartPointer<vtkIntArray>::New();
  bandPartNum->SetNumberOfComponents(1);
  bandPartNum->SetName("bandPartNum");
  
  vtkSmartPointer<vtkDoubleArray> bandPartNumAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandPartNumAVG->SetNumberOfComponents(1);
  bandPartNumAVG->SetName("bandPartNumAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandVol = vtkSmartPointer<vtkDoubleArray>::New();
  bandVol->SetNumberOfComponents(1);
  bandVol->SetName("bandVol");

  vtkSmartPointer<vtkDoubleArray> bandVolFraction = vtkSmartPointer<vtkDoubleArray>::New();
  bandVolFraction->SetNumberOfComponents(1);
  bandVolFraction->SetName("bandVolFraction");

  vtkSmartPointer<vtkDoubleArray> bandScherRate = vtkSmartPointer<vtkDoubleArray>::New();
  bandScherRate->SetNumberOfComponents(1);
  bandScherRate->SetName("bandScherRate");

  vtkSmartPointer<vtkDoubleArray> bandContactNumAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandContactNumAVG->SetNumberOfComponents(1);
  bandContactNumAVG->SetName("bandContactNumAVG");
  
  
  vtkSmartPointer<vtkDoubleArray> bandMu = vtkSmartPointer<vtkDoubleArray>::New();
  bandMu->SetNumberOfComponents(1);
  bandMu->SetName("bandMu");
  
  
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bandTMP = _bandRow->getBand(b);
    
    for (int z = 0; z<bandTMP->partNumb(); z++) {
      std::shared_ptr<particle> partTemp = bandTMP->getPart(z);
      if (not(partTemp->disabled())) {
        vtkIdType pid[1];
        pid[0] = spheresPos->InsertNextPoint(partTemp->c()[0], partTemp->c()[1], partTemp->c()[2]);
        radii->InsertNextValue(partTemp->rad());
        dist->InsertNextValue(partTemp->dist());
        height->InsertNextValue(partTemp->height());
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
        
        bandR->InsertNextValue(bandTMP->idR());
        bandZ->InsertNextValue(bandTMP->idZ());
        bandN->InsertNextValue(bandTMP->id());
        bandTau->InsertNextValue(bandTMP->tau());
        bandPress->InsertNextValue(bandTMP->press());
        bandLocalPress->InsertNextValue(bandTMP->localPress());
        bandLocSigmDev->InsertNextValue(bandTMP->dLocalAvg());
        bandOmega->InsertNextValue(bandTMP->omega());
        bandPartNum->InsertNextValue(bandTMP->partNumb());
        bandPartNumAVG->InsertNextValue(bandTMP->partNumb()/bandTMP->vol());
        bandVol->InsertNextValue(bandTMP->vol());
        bandVolFraction->InsertNextValue(bandTMP->volFraction());
        bandContactNumAVG->InsertNextValue(bandTMP->contactNumAVG());
        bandScherRate->InsertNextValue(bandTMP->scherRate());
        bandLocMu->InsertNextValue(bandTMP->muLocalAVG());
        //bandMu->InsertNextValue(bandTMP->tau()/bandTMP->press());
        
        spheresCells->InsertNextCell(1,pid);
      }
    }
    
    
    spheresUg->SetPoints(spheresPos);
    spheresUg->SetCells(VTK_VERTEX, spheresCells);
    spheresUg->GetPointData()->AddArray(radii);
    spheresUg->GetPointData()->AddArray(dist);
    spheresUg->GetPointData()->AddArray(height);
    spheresUg->GetPointData()->AddArray(spheresId);
    spheresUg->GetPointData()->AddArray(spheresType);
    spheresUg->GetPointData()->AddArray(spheresVelL);
    spheresUg->GetPointData()->AddArray(spheresVelA);
    spheresUg->GetPointData()->AddArray(vectorDr);
    spheresUg->GetPointData()->AddArray(vectorDz);
    spheresUg->GetPointData()->AddArray(vectorDf);
    spheresUg->GetPointData()->AddArray(bandR);
    spheresUg->GetPointData()->AddArray(bandZ);
    spheresUg->GetPointData()->AddArray(bandN);
    spheresUg->GetPointData()->AddArray(bandTau);
    spheresUg->GetPointData()->AddArray(bandPress);
    spheresUg->GetPointData()->AddArray(bandLocalPress);
    spheresUg->GetPointData()->AddArray(bandOmega);
    spheresUg->GetPointData()->AddArray(bandPartNum);
    spheresUg->GetPointData()->AddArray(bandPartNumAVG);
    spheresUg->GetPointData()->AddArray(bandVol);
    spheresUg->GetPointData()->AddArray(bandVolFraction);
    spheresUg->GetPointData()->AddArray(bandContactNumAVG);
    spheresUg->GetPointData()->AddArray(bandScherRate);
    spheresUg->GetPointData()->AddArray(bandLocSigmDev);
    spheresUg->GetPointData()->AddArray(bandLocMu);
    
  }
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetInput(spheresUg);
  writer->SetFileName("a.vtu");
  writer->Write();
};

void exportclass::gnuplotSchearRate() {
  
  
};
