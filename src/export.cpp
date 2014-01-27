/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RheometerAnalyze is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RheometerAnalyze.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "export.h"
#include <boost/foreach.hpp>
#include <cmath>

exportclass::exportclass(const std::shared_ptr<configopt> cfg, const std::shared_ptr <bandRow> bandAll, const std::vector<std::shared_ptr <forceRow>> forceAll) {
  _cfg = cfg;
  _bandRow = bandAll;
  _forceAll = forceAll;
};

void exportclass::VTK() {
  ofstream fileVTK;
  std::string _fileNameVTU;
  std::string _fileNameVTP;

  _fileNameVTU =  _cfg->FOutput();
  _fileNameVTP =  _cfg->FOutput();
  _fileNameVTU += "/output.vtu";
  _fileNameVTP += "/output.vtp";

  
  //Export Particles
  vtkSmartPointer<vtkPoints>  spheresPos = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkCellArray> spheresCells = vtkSmartPointer<vtkCellArray>::New();

  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("radii");

  vtkSmartPointer<vtkDoubleArray> mass = vtkSmartPointer<vtkDoubleArray>::New();
  mass->SetNumberOfComponents(1);
  mass->SetName("mass");

  vtkSmartPointer<vtkDoubleArray> density = vtkSmartPointer<vtkDoubleArray>::New();
  density->SetNumberOfComponents(1);
  density->SetName("density");
  
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
  
  vtkSmartPointer<vtkDoubleArray> posZyl = vtkSmartPointer<vtkDoubleArray>::New();
  posZyl->SetNumberOfComponents(3);
  posZyl->SetName("posZyl");
  
  vtkSmartPointer<vtkIntArray> bandR = vtkSmartPointer<vtkIntArray>::New();
  bandR->SetNumberOfComponents(1);
  bandR->SetName("bandR");
  
  vtkSmartPointer<vtkIntArray> bandZ = vtkSmartPointer<vtkIntArray>::New();
  bandZ->SetNumberOfComponents(1);
  bandZ->SetName("bandZ");
  
  vtkSmartPointer<vtkIntArray> bandF = vtkSmartPointer<vtkIntArray>::New();
  bandF->SetNumberOfComponents(1);
  bandF->SetName("bandF");
  
  vtkSmartPointer<vtkIntArray> bandN = vtkSmartPointer<vtkIntArray>::New();
  bandN->SetNumberOfComponents(1);
  bandN->SetName("bandN");
  
  vtkSmartPointer<vtkDoubleArray> bandTau = vtkSmartPointer<vtkDoubleArray>::New();
  bandTau->SetNumberOfComponents(1);
  bandTau->SetName("bandTau");
  
  vtkSmartPointer<vtkDoubleArray> bandPress = vtkSmartPointer<vtkDoubleArray>::New();
  bandPress->SetNumberOfComponents(1);
  bandPress->SetName("bandPress");
  
  vtkSmartPointer<vtkDoubleArray> bandTensor = vtkSmartPointer<vtkDoubleArray>::New();
  bandTensor->SetNumberOfComponents(9);
  bandTensor->SetName("bandTensor");
  
  vtkSmartPointer<vtkDoubleArray> partTensor = vtkSmartPointer<vtkDoubleArray>::New();
  partTensor->SetNumberOfComponents(9);
  partTensor->SetName("partTensor");
  
  vtkSmartPointer<vtkDoubleArray> bandTensorCap = vtkSmartPointer<vtkDoubleArray>::New();
  bandTensorCap->SetNumberOfComponents(9);
  bandTensorCap->SetName("bandTensorCap");
  
  vtkSmartPointer<vtkDoubleArray> partTensorCap = vtkSmartPointer<vtkDoubleArray>::New();
  partTensorCap->SetNumberOfComponents(9);
  partTensorCap->SetName("partTensorCap");
  
  vtkSmartPointer<vtkDoubleArray> bandMu = vtkSmartPointer<vtkDoubleArray>::New();
  bandMu->SetNumberOfComponents(1);
  bandMu->SetName("bandMu");
  
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
  
  vtkSmartPointer<vtkDoubleArray> bandVelLin = vtkSmartPointer<vtkDoubleArray>::New();
  bandVelLin->SetNumberOfComponents(3);
  bandVelLin->SetName("bandVelLin_rzf");
  
  vtkSmartPointer<vtkDoubleArray> bandWetContactsAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandWetContactsAVG->SetNumberOfComponents(1);
  bandWetContactsAVG->SetName("bandWetContactsAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandWetContactDistanceAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandWetContactDistanceAVG->SetNumberOfComponents(1);
  bandWetContactDistanceAVG->SetName("bandWetContactDistanceAVG");
  
  #ifdef ALGLIB
    vtkSmartPointer<vtkIntArray> bandShearBand = vtkSmartPointer<vtkIntArray>::New();
    bandShearBand->SetNumberOfComponents(1);
    bandShearBand->SetName("bandShearBand");
  #endif
  
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  
  if (_cfg->Vtk()==3 ) {                                // Force (interactions) to be exported
    vtkSmartPointer<vtkPoints> partPos = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> forceCells = vtkSmartPointer<vtkCellArray>::New();
    
    vtkSmartPointer<vtkFloatArray> force = vtkSmartPointer<vtkFloatArray>::New();
    force->SetNumberOfComponents(1);
    force->SetName("force");
    
    vtkSmartPointer<vtkIntArray> wet = vtkSmartPointer<vtkIntArray>::New();
    wet->SetNumberOfComponents(1);
    wet->SetName("wet");
    
    #ifdef ALGLIB
    vtkSmartPointer<vtkIntArray> shearband = vtkSmartPointer<vtkIntArray>::New();
    shearband->SetNumberOfComponents(1);
    shearband->SetName("shearband");
    #endif
    
    BOOST_FOREACH(std::shared_ptr <forceRow> fR,  _forceAll) {
      for(unsigned long long b=0; b<fR->arraySize(); b++) {
        
        partPos->InsertNextPoint(fR->getF(b)->pos1()(0), fR->getF(b)->pos1()(1), fR->getF(b)->pos1()(2));
        partPos->InsertNextPoint(fR->getF(b)->pos2()(0), fR->getF(b)->pos2()(1), fR->getF(b)->pos2()(2));
        
        vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
        line->GetPointIds()->SetId(0,(b*2));
        line->GetPointIds()->SetId(1,(b*2+1));
        forceCells->InsertNextCell(line);
        force->InsertNextValue(fR->getF(b)->val().norm());
        
        if (fR->getF(b)->volWater()>0) {
          wet->InsertNextValue(1);
        } else {
          wet->InsertNextValue(0);
        }
        
        #ifdef ALGLIB
        if (fR->getF(b)->part1()->shearBand() and fR->getF(b)->part2()->shearBand()) {
          shearband->InsertNextValue(2);
        } else if (fR->getF(b)->part1()->shearBand() or fR->getF(b)->part2()->shearBand()) {
          shearband->InsertNextValue(1);
        } else {
          shearband->InsertNextValue(0);
        }
        #endif
      }
    }
    
    vtkSmartPointer<vtkPolyData> fPd = vtkSmartPointer<vtkPolyData>::New();
    
    fPd->SetPoints(partPos);
    fPd->SetLines(forceCells);
    fPd->GetCellData()->AddArray(force);
    fPd->GetCellData()->AddArray(wet);
    #ifdef ALGLIB
    fPd->GetCellData()->AddArray(shearband);
    #endif
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    
    writer->SetDataModeToAscii();
    writer->SetInput(fPd);
    
    writer->SetFileName(_fileNameVTP.c_str());
    writer->Write();
    
  } else {
    for(unsigned int b=0; b<_bandRow->size(); b++) {
      std::shared_ptr<band> bandTMP = _bandRow->getBand(b);
      
      for (int z = 0; z<bandTMP->partNumb(); z++) {
        std::shared_ptr<particle> partTemp = bandTMP->getPart(z);
        if (not(partTemp->disabled())) {
          
          vtkIdType pid[1];
          if (_cfg->Vtk()==1) {       // Large VTK-file with all particles, snapshots
            pid[0] = spheresPos->InsertNextPoint(partTemp->c()[0], partTemp->c()[1], partTemp->c()[2]);
            radii->InsertNextValue(partTemp->rad());
            mass->InsertNextValue(partTemp->mass());
            density->InsertNextValue(partTemp->density());
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
            
            double posZ[3] = {partTemp->posZyl()(0), partTemp->posZyl()(1), partTemp->posZyl()(2)};
            posZyl->InsertNextTupleValue(posZ);
          } else if (_cfg->Vtk()==2 and (bandTMP->partNumb() > 0)) {    // Small VTK-file only with bands, averaged
            pid[0] = spheresPos->InsertNextPoint(bandTMP->midLinedR(), bandTMP->midLinedF(), bandTMP->midLinedZ());
            radii->InsertNextValue(bandTMP->radAvg());
            density->InsertNextValue(bandTMP->density());
            spheresType->InsertNextValue(bandTMP->type());
          }
          
          Eigen::Matrix3d tensorM = bandTMP->TensorAVG();
          
          double tensor[9] = {tensorM(0), tensorM(1), tensorM(2), 
                              tensorM(3), tensorM(4), tensorM(5), 
                              tensorM(6), tensorM(7), tensorM(8)};
          
          bandTensor->InsertNextTupleValue(tensor);
          
          tensorM = partTemp->stressTensorAVG();
          
          double  tensor2[9] = {tensorM(0), tensorM(1), tensorM(2), 
                              tensorM(3), tensorM(4), tensorM(5), 
                              tensorM(6), tensorM(7), tensorM(8)};
          
          partTensor->InsertNextTupleValue(tensor2);
          
          tensorM = bandTMP->TensorCapAVG();
          
          // Capillar tensor
          
          double tensor3[9] = {tensorM(0), tensorM(1), tensorM(2), 
                              tensorM(3), tensorM(4), tensorM(5), 
                              tensorM(6), tensorM(7), tensorM(8)};
          
          bandTensorCap->InsertNextTupleValue(tensor3);
          
          tensorM = partTemp->stressTensorCapAVG();
          
          double  tensor4[9] = {tensorM(0), tensorM(1), tensorM(2), 
                              tensorM(3), tensorM(4), tensorM(5), 
                              tensorM(6), tensorM(7), tensorM(8)};
          
          partTensorCap->InsertNextTupleValue(tensor4);
          
          bandR->InsertNextValue(bandTMP->idR());
          bandZ->InsertNextValue(bandTMP->idZ());
          bandF->InsertNextValue(bandTMP->idF());
          bandN->InsertNextValue(bandTMP->id());
          bandTau->InsertNextValue(bandTMP->tau());
          bandPress->InsertNextValue(bandTMP->press());
          bandOmega->InsertNextValue(bandTMP->omega());
          bandPartNum->InsertNextValue(bandTMP->partNumb());
          bandPartNumAVG->InsertNextValue(bandTMP->partNumb()/bandTMP->vol());
          bandVol->InsertNextValue(bandTMP->vol());
          bandVolFraction->InsertNextValue(bandTMP->volFraction());
          bandContactNumAVG->InsertNextValue(bandTMP->contactNumAVG());
          bandScherRate->InsertNextValue(bandTMP->scherRate());
          bandMu->InsertNextValue(bandTMP->muAVG());
          bandWetContactDistanceAVG->InsertNextValue(bandTMP->wetContactDistanceAVG());
          bandWetContactsAVG->InsertNextValue(bandTMP->wetContactsAVG());
          
          double VelLin[3] = {bandTMP->vZyl()[0], bandTMP->vZyl()[1], bandTMP->vZyl()[2]};
          bandVelLin->InsertNextTupleValue(VelLin);
          #ifdef ALGLIB
            if (bandTMP->shearBand()) {
              bandShearBand->InsertNextValue(1);
            } else {
              bandShearBand->InsertNextValue(0);
            }
          #endif
          
          spheresCells->InsertNextCell(1,pid);
          if (_cfg->Vtk()==2) break;
        }
      }
      
      
      
      spheresUg->SetPoints(spheresPos);
      spheresUg->SetCells(VTK_VERTEX, spheresCells);
      spheresUg->GetPointData()->AddArray(radii);
      spheresUg->GetPointData()->AddArray(density);
      
      if (_cfg->Vtk()==1) {
        // Only for particles, _cfg->Vtk()==1
        spheresUg->GetPointData()->AddArray(mass);
        spheresUg->GetPointData()->AddArray(spheresId);
        spheresUg->GetPointData()->AddArray(spheresVelL);
        spheresUg->GetPointData()->AddArray(spheresVelA);
        spheresUg->GetPointData()->AddArray(vectorDr);
        spheresUg->GetPointData()->AddArray(vectorDz);
        spheresUg->GetPointData()->AddArray(vectorDf);
        spheresUg->GetPointData()->AddArray(posZyl);
      }
      
      spheresUg->GetPointData()->AddArray(spheresType);
      spheresUg->GetPointData()->AddArray(bandR);
      spheresUg->GetPointData()->AddArray(bandZ);
      spheresUg->GetPointData()->AddArray(bandF);
      spheresUg->GetPointData()->AddArray(bandN);
      spheresUg->GetPointData()->AddArray(bandTau);
      spheresUg->GetPointData()->AddArray(bandPress);
      spheresUg->GetPointData()->AddArray(bandOmega);
      spheresUg->GetPointData()->AddArray(bandPartNum);
      spheresUg->GetPointData()->AddArray(bandPartNumAVG);
      spheresUg->GetPointData()->AddArray(bandVol);
      spheresUg->GetPointData()->AddArray(bandVolFraction);
      spheresUg->GetPointData()->AddArray(bandContactNumAVG);
      spheresUg->GetPointData()->AddArray(bandScherRate);
      spheresUg->GetPointData()->AddArray(bandMu);
      spheresUg->GetPointData()->AddArray(bandVelLin);
      spheresUg->GetPointData()->AddArray(bandTensor);
      spheresUg->GetPointData()->AddArray(partTensor);
      spheresUg->GetPointData()->AddArray(bandTensorCap);
      spheresUg->GetPointData()->AddArray(partTensorCap);
      spheresUg->GetPointData()->AddArray(bandWetContactsAVG);
      spheresUg->GetPointData()->AddArray(bandWetContactDistanceAVG);
      #ifdef ALGLIB
      spheresUg->GetPointData()->AddArray(bandShearBand);
      #endif
    }
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetDataModeToAscii();
    writer->SetInput(spheresUg);
    
    writer->SetFileName(_fileNameVTU.c_str());
    writer->Write();
  }
};

void exportclass::gnuplotSchearRate() {  
  
  std::string _fileNameGS, _fileNameG;

  _fileNameGS =  _cfg->FOutput();
  _fileNameGS += "/shearband";
  
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/gnuplot_daten";
  
  ofstream myfile4 (_fileNameGS.c_str());
  #ifdef ALGLIB
  if (myfile4.is_open()) {
    myfile4 << "RZ\tW\tH" << std::endl;
    for(unsigned int b=0; b<_bandRow->shearBandSize(); b++) {
      const Eigen::Vector3d tmpVal = _bandRow->shearBand(b);
      myfile4 << tmpVal[0] << "\t"<< tmpVal[1] << "\t"<< tmpVal[2] << std::endl;
    }
  }
  #endif
  
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_id\t002_r\t003_z\t004_rPos\t005_zPos\t006_H\t007_W\t008_PartN\t";
  myfileG << "009_vol\t010_volFract\t011_contactNum\t012_vZyldR\t013_vZyldZ\t014_vZyldF\t";
  myfileG << "015_ShearRate\t016_I\t017_Eta\t018_omega\t019_press\t020_tau\t021_mu\t";
  myfileG << "022_strTensRR\t023_strTensRZ\t024_strTensRF\t";
  myfileG << "025_strTensZR\t026_strTensZZ\t027_strTensZF\t";
  myfileG << "028_strTensFR\t029_strTensFZ\t030_strTensFF\t";
  myfileG << "031_ShearBand\t032_WetContactsAVG\t033_WetContactsDistAVG\t";
  myfileG << "034_f\t035_fPos\t";
  myfileG << "036_strTensCapRR\t037_strTensCapRZ\t038_strTensCapRF\t";
  myfileG << "039_strTensCapZR\t040_strTensCapZZ\t041_strTensCapZF\t";
  myfileG << "042_strTensCapFR\t043_strTensCapFZ\t044_strTensCapFF\t \n";
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bT = _bandRow->getBand(b);
    myfileG << bT->id() << "\t";           // 001_id
    myfileG << bT->idR() << "\t";          // 002_r
    myfileG << bT->idZ() << "\t";          // 003_z
    myfileG << bT->midLinedR() << "\t";    // 004_rPos
    myfileG << bT->midLinedZ() << "\t";    // 005_zPos
    myfileG << bT->dZ() << "\t";           // 006_H
    myfileG << bT->dR() << "\t";           // 007_W
    myfileG << bT->partNumb() << "\t";     // 008_PartN
    myfileG << bT->vol() << "\t";          // 009_vol
    myfileG << bT->volFraction() << "\t";  // 010_volFract
    myfileG << bT->contactNumAVG() << "\t";// 011_contactNum
    myfileG << bT->vZyl()[0]<< "\t";       // 012_vZyldR
    myfileG << bT->vZyl()[1]<< "\t";       // 013_vZyldZ
    myfileG << bT->vZyl()[2]<< "\t";       // 014_vZyldF
    myfileG << bT->scherRate()<< "\t";     // 015_ShearRate
    myfileG << bT->I()<< "\t";             // 016_I
    myfileG << bT->eta()<< "\t";           // 017_Eta
    myfileG << bT->omega()<< "\t";         // 018_omega
    myfileG << bT->press()<< "\t";         // 019_press
    myfileG << bT->tau()<< "\t";           // 020_tau
    myfileG << bT->muAVG()<< "\t";         // 021_mu
    myfileG << bT->TensorAVG()(0)<< "\t";  // 022_strTensRR
    myfileG << bT->TensorAVG()(1)<< "\t";  // 023_strTensRZ
    myfileG << bT->TensorAVG()(2)<< "\t";  // 024_strTensRF
    myfileG << bT->TensorAVG()(3)<< "\t";  // 025_strTensZR
    myfileG << bT->TensorAVG()(4)<< "\t";  // 026_strTensZZ
    myfileG << bT->TensorAVG()(5)<< "\t";  // 027_strTensZF
    myfileG << bT->TensorAVG()(6)<< "\t";  // 028_strTensFR
    myfileG << bT->TensorAVG()(7)<< "\t";  // 029_strTensFZ
    myfileG << bT->TensorAVG()(8)<< "\t";  // 030_strTensFF
    myfileG << bT->shearBand()<< "\t";     // 031_ShearBand - True (1) or False (0)
    myfileG << bT->wetContactsAVG()<< "\t";// 032_WetContactsAVG
    myfileG << bT->wetContactDistanceAVG()<< "\t";// 033_WetContactsDistAVG
    myfileG << bT->idF() << "\t";             // 034_f
    myfileG << bT->midLinedF() << "\t";       // 035_fPos
    myfileG << bT->TensorCapAVG()(0)<< "\t";  // 036_strTensCapRR
    myfileG << bT->TensorCapAVG()(1)<< "\t";  // 037_strTensCapRZ
    myfileG << bT->TensorCapAVG()(2)<< "\t";  // 038_strTensCapRF
    myfileG << bT->TensorCapAVG()(3)<< "\t";  // 039_strTensCapZR
    myfileG << bT->TensorCapAVG()(4)<< "\t";  // 040_strTensCapZZ
    myfileG << bT->TensorCapAVG()(5)<< "\t";  // 041_strTensCapZF
    myfileG << bT->TensorCapAVG()(6)<< "\t";  // 042_strTensCapFR
    myfileG << bT->TensorCapAVG()(7)<< "\t";  // 043_strTensCapFZ
    myfileG << bT->TensorCapAVG()(8)<< "\t";  // 044_strTensCapFF
    myfileG << " \n"; 
  }
        
  myfileG.close();
};

void exportclass::gnuplotContactAnalyze(int bins) {  
  // Calculates, how many contacts are having the corresponding distance Delta
  std::string _fileNameG;
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/contacts";
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_id\t002_minDelta\t003_maxDelta\t004_ContNumber\t005_ContNumberAVG\t006_ForceAVG\n";
  
  std::vector <std::shared_ptr<force> >  deltas;
  std::vector <long long int>            deltasBin(bins);
  std::vector <double>                   forcesBin(bins);
  
  for(int x = 0; x < bins; ++x) {
    deltasBin[x] = 0;
    forcesBin[x] = 0;
  }
  
  double minDelta = 0.0;
  double maxDelta = 0.0;
    
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<force> > forces = snapshotCur->forces();
    
    BOOST_FOREACH(std::shared_ptr<force> f, forces) {
      long long pid1T = f->pid1();
      long long pid2T = f->pid2();
      
      if ((_cfg->tF()>=0) and (f->part1()->type() != _cfg->tF()) and (f->part2()->type() == _cfg->tF())) {
        pid1T = -1;
      }
      if ((_cfg->tF()>=0) and (f->part2()->type() != _cfg->tF()) and (f->part1()->type() == _cfg->tF())) {
        pid2T = -1;
      }
      
      if ((pid1T>=0) and (pid2T>=0)) {
        deltas.push_back(f);
        if  (minDelta == maxDelta and minDelta == 0 ) {
          minDelta = f->deltaN();
          maxDelta = f->deltaN();
        }
        minDelta = std::min(minDelta, f->deltaN());
        maxDelta = std::max(maxDelta, f->deltaN());
      }
    }
  }
  
  double DDelta = (maxDelta - minDelta)/bins;
  BOOST_FOREACH(std::shared_ptr<force> d, deltas) {
    deltasBin[int(floor((d->deltaN()-minDelta)/DDelta))] += 1;
    forcesBin[int(floor((d->deltaN()-minDelta)/DDelta))] += d->val().norm();
  }
  
  for(unsigned int x = 0; x < deltasBin.size(); ++x) {
    double forceTmp = 0.0;
    if (deltasBin[x]!=0) {
      forceTmp = forcesBin[x]/deltasBin[x];
    }
    myfileG << x << " " <<minDelta + DDelta*x  << " " <<minDelta + DDelta*(x+1) 
            << " " <<  deltasBin[x] << " " <<  deltasBin[x]/snapshots->size() 
            << " "<< forceTmp <<  "\n";
  }
  
  myfileG.close();
};

void exportclass::gnuplotContactWet() {  
  // Calculates, how many contacts are having the corresponding distance Delta
  std::string _fileNameG;
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/contactsWet";
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_iter\t002_time\t003_Contacts\t004_WetContacts\t005_LiquidVol\n";
  
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<force> > forces = snapshotCur->forces();
    unsigned long long int numContacts = 0;
    unsigned long long int numWetContacts = 0;
    double liqVol = 0;
    
    BOOST_FOREACH(std::shared_ptr<force> f, forces) {
      long long pid1T = f->pid1();
      long long pid2T = f->pid2();
      
      if ((_cfg->tF()>=0) and (f->part1()->type() != _cfg->tF()) and (f->part2()->type() == _cfg->tF())) {
        pid1T = -1;
      }
      if ((_cfg->tF()>=0) and (f->part2()->type() != _cfg->tF()) and (f->part1()->type() == _cfg->tF())) {
        pid2T = -1;
      }
      
      if ((pid1T>=0) and (pid2T>=0)) {
        numContacts++;
        if (f->volWater()>0) {
          numWetContacts++;
          liqVol+=f->volWater();
        }
      }
    }
    
    myfileG << snapshotCur->timeStep() << " " << snapshotCur->timeStep()*_cfg->dT()  << " " << numContacts << " " <<  numWetContacts << " " <<  liqVol  << "\n";
    
  }
  myfileG.close();
};


bool sortContactFollow(std::shared_ptr<contactFollow> i, std::shared_ptr<contactFollow> j) {
  
  if (i->P1_id() < j->P1_id()) {
    return true;
  } else {
    if (i->P2_id() < j->P2_id() and 
        i->P1_id() == j->P1_id()) {
      return true;
    } else {
      if (i->timeStep() < j->timeStep() and 
          i->P2_id() == j->P2_id() and 
          i->P1_id() == j->P1_id()) {
        return true;
      } else {
        return false;
      }
    }
  }
}

void exportclass::gnuplotContactFollow() {
  std::string _fileNameG;
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/followContacts";
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_Pid1\t002_Pid2\t003_time\t004_delta\t005_force\t006_volWater\t007_distCurr\t008_distCrit\n";
  
  
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  std::vector <std::shared_ptr<contactFollow> > cntFolw;
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<force> > forces = snapshotCur->forces();
    
    BOOST_FOREACH(std::shared_ptr<force> f, forces) {
      std::shared_ptr<contactFollow> tmpCntFolw ( new contactFollow (f, snapshotCur));
      cntFolw.push_back(tmpCntFolw);
    }
  }
  
  std::sort (cntFolw.begin(), cntFolw.end(), sortContactFollow);
  
  BOOST_FOREACH(std::shared_ptr<contactFollow> tmpCntFolw, cntFolw) {
     myfileG << tmpCntFolw->P1_id() << " " << tmpCntFolw->P2_id() << " " 
             << tmpCntFolw->timeStep()*_cfg->dT() << " " 
             << tmpCntFolw->deltaN() << " " 
             << tmpCntFolw->f().norm() << " "
             << tmpCntFolw->volWater() << " "
             << tmpCntFolw->distCurr() << " "
             << tmpCntFolw->distCrit() << " " << std::endl;
   }
  
  myfileG.close();
};



void exportclass::Utwente()  {
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  //unsigned int leadingZerosLen = static_cast <unsigned int> (log10 (snapshots->size()) + 1);
  unsigned int leadingZerosLen = 4;
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    
    std::stringstream ss;
    ss << setw(leadingZerosLen) << setfill('0') << i;
    
    double timeTmp = snapshotCur->timeStep()*_cfg->dT();
    //================================================
    
    std::string _fileNameC3d;
    _fileNameC3d = _cfg->FOutput();
    _fileNameC3d += "/c3d.";
    _fileNameC3d += ss.str();
    
    ofstream C3d (_fileNameC3d.c_str());
    std::vector <std::shared_ptr<particle> > particles = snapshotCur->particles();
    C3d << particles.size() << " " << timeTmp
        << " " << -_cfg->Dout()/2.0 << " "<<-_cfg->Dout()/2.0 << " 0.0 "
        << _cfg->Dout()/2.0 << " " << _cfg->Dout()/2.0 << " " << _cfg->H()
        << std::endl;
    
    BOOST_FOREACH(std::shared_ptr<particle> p, particles) {
      C3d << p->c()(0) << " " << p->c()(1) << " "<< p->c()(2) << " " 
          << p->v()(0) << " " << p->v()(1) << " "<< p->v()(2) << " "
          << p->rad() << " 0 " << std::endl;
    };
    
    C3d << "\t"<< std::endl;
    C3d.close();
    
    //================================================
    std::string _fileNameFstat;
    _fileNameFstat = _cfg->FOutput();
    _fileNameFstat += "/fstat.";
    _fileNameFstat += ss.str();
    
    ofstream Fstat (_fileNameFstat.c_str());
    std::vector <std::shared_ptr<force> > forces = snapshotCur->forces();
    Fstat<<"# " << timeTmp << " volumenfraction" << std::endl;
    Fstat<<"# Time " << timeTmp << std::endl;
    Fstat<<"#  " << std::endl;
    
    
    BOOST_FOREACH(std::shared_ptr<force> f, forces) {
      long long pid1T = f->pid1();
      long long pid2T = f->pid2();
      
      if ((_cfg->tF()>=0) and (f->part1()->type() != _cfg->tF()) and (f->part2()->type() == _cfg->tF())) {
        pid1T = -1;
      }
      if ((_cfg->tF()>=0) and (f->part2()->type() != _cfg->tF()) and (f->part1()->type() == _cfg->tF())) {
        pid2T = -1;
      }
      
      if (pid1T>=0) {
        Fstat << timeTmp << " " << pid1T << " " << pid2T << " "
              << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
              << f->deltaN() << " "
              << "0.0" << " "                                                 //Theta, should be fixed.
              << (f->forceN()).norm() << " " <<  (f->forceT()).norm() << " "  //Normal and tangential components of the force
              << f->nVec()(0) << " " << f->nVec()(1) << " " << f->nVec()(2) << " " 
              << f->tVec()(0) << " " << f->tVec()(1) << " " << f->tVec()(2) << " " 
              << std::endl;
      } else {
        Fstat << timeTmp << " " << pid2T << " " << pid1T << " "
              << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
              << f->deltaN() << " "
              << "0.0" << " "                                                 //Theta, should be fixed.
              << (f->forceN()).norm() << " " <<  (f->forceT()).norm() << " "  //Normal and tangential components of the force
              << -f->nVec()(0) << " " << -f->nVec()(1) << " " << -f->nVec()(2) << " " 
              << -f->tVec()(0) << " " << -f->tVec()(1) << " " << -f->tVec()(2) << " " 
              << std::endl;
      }
      if ((pid1T>=0) and (pid2T>=0)) {
        Fstat << timeTmp << " " << pid2T << " " << pid1T << " "
              << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
              << f->deltaN() << " " 
              << "0.0" << " "                                                 //Theta, should be fixed.
              << (f->forceN()).norm() << " " <<  (f->forceT()).norm() << " "  //Normal and tangential components of the force
              << -f->nVec()(0) << " " << -f->nVec()(1) << " " << -f->nVec()(2) << " " 
              << -f->tVec()(0) << " " << -f->tVec()(1) << " " << -f->tVec()(2) << " " 
              << std::endl;
      }
    }
    
    Fstat.close();
  }
}

void exportclass::intOri() {  
  // Exports orientations of interactions
  std::string _fileName;
  _fileName  =  _cfg->FOutput();
  _fileName  +=  "/interOri";
  
  ofstream myfile (_fileName.c_str());
  myfile << "#001_id\t002_r\t003_z\t004_rPos\t005_zPos\t006_Theta\t007_Psi\t008_normIterOri\t009_normIterOriN\t";
  myfile << "#010_capiIterOri\t#011_capiIterOriN\t";
  myfile << "#012_SphX\t#013_SphY#014_SphZ\t\n";
  myfile << "#015_SphXnormContOriN[Phi]\t#016_SphYnormContOriN[Rho]#017_SphZnormContOriN[Z]\t\n";
  myfile << "#018_SphXcapiContOriN[Phi]\t#019_SphYcapiContOriN[Rho]#020_SphZcapiContOriN[Z]\t\n";
  
  
  unsigned long numbLine=0;
  const double dAngle = M_PI/_cfg->intOri();
  const double d2Ang  = dAngle/2.0;
  
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bT = _bandRow->getBand(b);
    InteractionsMatrixD  normContOri = bT->normContOri();
    InteractionsMatrixD  normContOriN = normContOri;
    InteractionsMatrixD  capiContOri = bT->capiContOri();
    InteractionsMatrixD  capiContOriN = capiContOri;
    
    
    if (normContOriN.norm()>0.0){
      normContOriN.normalize();                     // Normalized normal contact number in every slot
    }

    if (capiContOriN.norm()>0.0){
      capiContOriN.normalize();                     // Normalized capillary contact number in every slot
    }
    
    for (unsigned short ThetaI=0; ThetaI < _cfg->intOri()*2; ThetaI++) {
      for (unsigned short PsiI=0; PsiI   < _cfg->intOri(); PsiI++) {
        Eigen::Vector3d SphXYZ = sph_to_cart(Eigen::Vector3d(ThetaI*dAngle + d2Ang, PsiI*dAngle + d2Ang, 1.0));
        const double ThetaId = ThetaI*dAngle + d2Ang;
        const double PsiId   = PsiI*dAngle + d2Ang;
        
        Eigen::Vector3d SphXYZ_normContOriN = sph_to_cart(Eigen::Vector3d(ThetaId, PsiId, normContOriN(ThetaI, PsiI)));
        Eigen::Vector3d SphXYZ_capiContOriN = sph_to_cart(Eigen::Vector3d(ThetaId, PsiId, capiContOriN(ThetaI, PsiI)));
        Eigen::Vector3d SphXYZ_normContOri  = sph_to_cart(Eigen::Vector3d(ThetaId, PsiId, normContOri(ThetaI, PsiI)));
        Eigen::Vector3d SphXYZ_capiContOri  = sph_to_cart(Eigen::Vector3d(ThetaId, PsiId, capiContOri(ThetaI, PsiI)));
        
        myfile << bT->id() << "\t";                   // 001_id
        myfile << bT->idR() << "\t";                  // 002_r
        myfile << bT->idZ() << "\t";                  // 003_z
        myfile << bT->midLinedR() << "\t";            // 004_rPos
        myfile << bT->midLinedZ() << "\t";            // 005_zPos
        myfile << ThetaI*dAngle + d2Ang << "\t";      // 006_Theta
        myfile << PsiI*dAngle + d2Ang << "\t";        // 007_Psi
        myfile << normContOri(ThetaI, PsiI) << "\t";  // 008_normIterOri
        myfile << normContOriN(ThetaI, PsiI) << "\t"; // 009_normIterOriN
        myfile << capiContOri(ThetaI, PsiI) << "\t";  // 010_capiIterOri
        myfile << capiContOriN(ThetaI, PsiI) << "\t"; // 011_capiIterOriN
        myfile << SphXYZ(0) << "\t";                  // 012_SphX
        myfile << SphXYZ(1) << "\t";                  // 013_SphY
        myfile << SphXYZ(2) << "\t";                  // 014_SphZ
        myfile << SphXYZ_normContOriN(0) << "\t";     // 015_SphXnormContOriN
        myfile << SphXYZ_normContOriN(1) << "\t";     // 016_SphYnormContOriN
        myfile << SphXYZ_normContOriN(2) << "\t";     // 017_SphZnormContOriN
        myfile << SphXYZ_capiContOriN(0) << "\t";     // 018_SphXcapiContOriN
        myfile << SphXYZ_capiContOriN(1) << "\t";     // 019_SphYcapiContOriN
        myfile << SphXYZ_capiContOriN(2) << "\t";     // 020_SphZcapiContOriN
        myfile << SphXYZ_normContOri(0)  << "\t";     // 021_SphXnormContOri
        myfile << SphXYZ_normContOri(1)  << "\t";     // 022_SphYnormContOri
        myfile << SphXYZ_normContOri(2)  << "\t";     // 023_SphZnormContOri
        myfile << SphXYZ_capiContOri(0)  << "\t";     // 024_SphXcapiContOri
        myfile << SphXYZ_capiContOri(1)  << "\t";     // 025_SphYcapiContOri
        myfile << SphXYZ_capiContOri(2)  << "\t";     // 026_SphZcapiContOri
        
        myfile << "\n";    // 
        numbLine++;
      }
    }
  }
  myfile.close();
  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  _fileName  =  _cfg->FOutput();
  _fileName  +=  "/interOri2D";
  
  ofstream myfile2 (_fileName.c_str());
  
  myfile2 << "#001_id\t002_r\t003_z\t004_rPos\t005_zPos\t006_Theta\t007_Psi\t008_normIterOri\t009_normIterOriN\t";
  myfile2 << "#010_capiIterOri\t#011_capiIterOriN\t";
  myfile2 << "#012_SphX\t#013_SphY#014_\t\n";
  myfile2 << "#015_SphXnormContOriN[Phi]\t#016_SphYnormContOriN[Rho]#017_\t\n";
  myfile2 << "#018_SphXcapiContOriN[Phi]\t#019_SphYcapiContOriN[Rho]#020_\t\n";
  
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bT = _bandRow->getBand(b);
    InteractionsMatrixD  normContOri    = bT->normContOri().colwise().sum();
    InteractionsMatrixD  normContOriN   = normContOri.colwise().sum();
    InteractionsMatrixD  capiContOri    = bT->capiContOri().colwise().sum();
    InteractionsMatrixD  capiContOriN   = capiContOri.colwise().sum();
    
    
    if (normContOriN.norm()>0.0){
      normContOriN.normalize();                     // Normalized normal contact number in every slot
    }

    if (capiContOriN.norm()>0.0){
      capiContOriN.normalize();                     // Normalized capillary contact number in every slot
    }
    
    unsigned short ThetaI = M_PI/2.0;
    
    for (unsigned short PsiI=0; PsiI   < _cfg->intOri(); PsiI++) {
      const double ThetaId = ThetaI*dAngle + d2Ang;
      const double PsiId   = PsiI*dAngle + d2Ang - M_PI/2.0;
      Eigen::Vector2d SphXYZ = Eigen::Vector2d(cos(PsiId), sin(PsiId));
      
      Eigen::Vector2d SphXYZ_normContOriN = Eigen::Vector2d(cos(PsiId)*normContOriN(PsiI,  0), sin(PsiId)*normContOriN(PsiI,  0));
      Eigen::Vector2d SphXYZ_capiContOriN = Eigen::Vector2d(cos(PsiId)*capiContOriN(PsiI,  0), sin(PsiId)*capiContOriN(PsiI,  0));
      Eigen::Vector2d SphXYZ_normContOri  = Eigen::Vector2d(cos(PsiId)*normContOri (PsiI,  0), sin(PsiId)*normContOri (PsiI,  0));
      Eigen::Vector2d SphXYZ_capiContOri  = Eigen::Vector2d(cos(PsiId)*capiContOri (PsiI,  0), sin(PsiId)*capiContOri (PsiI,  0));
      
      myfile2 << bT->id() << "\t";                   // 001_id
      myfile2 << bT->idR() << "\t";                  // 002_r
      myfile2 << bT->idZ() << "\t";                  // 003_z
      myfile2 << bT->midLinedR() << "\t";            // 004_rPos
      myfile2 << bT->midLinedZ() << "\t";            // 005_zPos
      myfile2 << ThetaI*dAngle + d2Ang << "\t";      // 006_Theta
      myfile2 << PsiId               << "\t";        // 007_Psi
      myfile2 << normContOri(ThetaI, PsiI) << "\t";  // 008_normIterOri
      myfile2 << normContOriN(ThetaI, PsiI) << "\t"; // 009_normIterOriN
      myfile2 << capiContOri(ThetaI, PsiI) << "\t";  // 010_capiIterOri
      myfile2 << capiContOriN(ThetaI, PsiI) << "\t"; // 011_capiIterOriN
      myfile2 << SphXYZ(0) << "\t";                  // 012_SphX
      myfile2 << SphXYZ(1) << "\t";                  // 013_SphY
      myfile2              << "0.0\t";               // 014_SphZ
      myfile2 << SphXYZ_normContOriN(0) << "\t";     // 015_SphXnormContOriN
      myfile2 << SphXYZ_normContOriN(1) << "\t";     // 016_SphYnormContOriN
      myfile2                           << "0.0\t";  // 017_SphZnormContOriN
      myfile2 << SphXYZ_capiContOriN(0) << "\t";     // 018_SphXcapiContOriN
      myfile2 << SphXYZ_capiContOriN(1) << "\t";     // 019_SphYcapiContOriN
      myfile2                           << "0.0\t";  // 020_SphZcapiContOriN
      myfile2 << SphXYZ_normContOri(0)  << "\t";     // 021_SphXnormContOri
      myfile2 << SphXYZ_normContOri(1)  << "\t";     // 022_SphYnormContOri
      myfile2                           << "0.0\t";  // 023_SphZnormContOri
      myfile2 << SphXYZ_capiContOri(0)  << "\t";     // 024_SphXcapiContOri
      myfile2 << SphXYZ_capiContOri(1)  << "\t";     // 025_SphYcapiContOri
      myfile2                           << "0.0\t";  // 026_SphZcapiContOri
      
      myfile2 << "\n";    // 
      numbLine++;
    }
  }
  myfile2.close();
};

void exportclass::torque() {
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  std::string _fileName  =  _cfg->FOutput();
  _fileName  +=  "/torque";
  ofstream myfile2 (_fileName.c_str());
  
  myfile2 << "#001_Id\t002_time\t003_torque\n";
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    myfile2                                                             << i << "\t";  // 001_snapId
    myfile2                            << snapshotCur->timeStep()*_cfg->dT() << "\t";  // 002_time
    myfile2 << snapshotCur->torque(_cfg->get_o(), _cfg->get_c(), _cfg->tR()) << "\n";  // 003_torque
  }
  myfile2.close();
}
