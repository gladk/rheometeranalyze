/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014, 2015 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014, 2015 Anton Gladky <gladky.anton@gmail.com>

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
#include <boost/unordered_map.hpp>

exportclass::exportclass(const std::shared_ptr<configopt> cfg, const std::shared_ptr <bandRow> bandAll, const std::vector<std::shared_ptr <forceRow>> forceAll) {
  _cfg = cfg;
  _bandRow = bandAll;
  _forceAll = forceAll;
  _noForces = false;
};

exportclass::exportclass(const std::shared_ptr<configopt> cfg, const std::shared_ptr <bandRow> bandAll) {
  _cfg = cfg;
  _bandRow = bandAll;
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

  // Parameters per particle=====================================================
  vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
  radii->SetNumberOfComponents(1);
  radii->SetName("p_Radii");

  vtkSmartPointer<vtkDoubleArray> mass = vtkSmartPointer<vtkDoubleArray>::New();
  mass->SetNumberOfComponents(1);
  mass->SetName("p_Mass");

  vtkSmartPointer<vtkDoubleArray> density = vtkSmartPointer<vtkDoubleArray>::New();
  density->SetNumberOfComponents(1);
  density->SetName("p_Density");
  
  vtkSmartPointer<vtkIntArray> spheresId = vtkSmartPointer<vtkIntArray>::New();
  spheresId->SetNumberOfComponents(1);
  spheresId->SetName("p_Id");

  vtkSmartPointer<vtkIntArray> spheresType = vtkSmartPointer<vtkIntArray>::New();
  spheresType->SetNumberOfComponents(1);
  spheresType->SetName("p_Type");
  
  vtkSmartPointer<vtkIntArray> sphereSnapshot = vtkSmartPointer<vtkIntArray>::New();
  sphereSnapshot->SetNumberOfComponents(1);
  sphereSnapshot->SetName("p_Snapshot");
  
  vtkSmartPointer<vtkIntArray> sphereHighStress = vtkSmartPointer<vtkIntArray>::New();
  sphereHighStress->SetNumberOfComponents(1);
  sphereHighStress->SetName("p_Highstress");
  
  vtkSmartPointer<vtkDoubleArray> spheresVelL = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelL->SetNumberOfComponents(3);
  spheresVelL->SetName("p_Velocity_lin");

  vtkSmartPointer<vtkDoubleArray> spheresVelA = vtkSmartPointer<vtkDoubleArray>::New();
  spheresVelA->SetNumberOfComponents(3);
  spheresVelA->SetName("p_Velocity_ang");
  
  vtkSmartPointer<vtkDoubleArray> vectorDr = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDr->SetNumberOfComponents(3);
  vectorDr->SetName("p_Vector_dr");
  
  vtkSmartPointer<vtkDoubleArray> vectorDz = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDz->SetNumberOfComponents(3);
  vectorDz->SetName("p_Vector_dz");
  
  vtkSmartPointer<vtkDoubleArray> vectorDf = vtkSmartPointer<vtkDoubleArray>::New();
  vectorDf->SetNumberOfComponents(3);
  vectorDf->SetName("p_Vector_df");
  
  vtkSmartPointer<vtkDoubleArray> posZyl = vtkSmartPointer<vtkDoubleArray>::New();
  posZyl->SetNumberOfComponents(3);
  posZyl->SetName("p_PosZyl");
  
  vtkSmartPointer<vtkDoubleArray> bandTensorCap = vtkSmartPointer<vtkDoubleArray>::New();
  bandTensorCap->SetNumberOfComponents(9);
  bandTensorCap->SetName("b_TensorCap");
  
  vtkSmartPointer<vtkDoubleArray> partTensorCap = vtkSmartPointer<vtkDoubleArray>::New();
  partTensorCap->SetNumberOfComponents(9);
  partTensorCap->SetName("p_TensorCap");
  
  vtkSmartPointer<vtkDoubleArray> stress3 = vtkSmartPointer<vtkDoubleArray>::New();
  stress3->SetNumberOfComponents(1);
  stress3->SetName("p_Stress3");
  
  vtkSmartPointer<vtkDoubleArray> stress1 = vtkSmartPointer<vtkDoubleArray>::New();
  stress1->SetNumberOfComponents(1);
  stress1->SetName("p_Stress1");
  
  vtkSmartPointer<vtkIntArray> contactNum = vtkSmartPointer<vtkIntArray>::New();
  contactNum->SetNumberOfComponents(1);
  contactNum->SetName("p_ContactNum");
  
  vtkSmartPointer<vtkIntArray> wetContactNum = vtkSmartPointer<vtkIntArray>::New();
  wetContactNum->SetNumberOfComponents(1);
  wetContactNum->SetName("p_WetContactNum");
  
  vtkSmartPointer<vtkDoubleArray> volWater = vtkSmartPointer<vtkDoubleArray>::New();
  volWater->SetNumberOfComponents(1);
  volWater->SetName("p_VolWater");
  
  // Parameters per band=====================================================
  
  vtkSmartPointer<vtkIntArray> bandR = vtkSmartPointer<vtkIntArray>::New();
  bandR->SetNumberOfComponents(1);
  bandR->SetName("b_R");
  
  vtkSmartPointer<vtkIntArray> bandZ = vtkSmartPointer<vtkIntArray>::New();
  bandZ->SetNumberOfComponents(1);
  bandZ->SetName("b_Z");
  
  vtkSmartPointer<vtkIntArray> bandF = vtkSmartPointer<vtkIntArray>::New();
  bandF->SetNumberOfComponents(1);
  bandF->SetName("b_F");
  
  vtkSmartPointer<vtkIntArray> bandN = vtkSmartPointer<vtkIntArray>::New();
  bandN->SetNumberOfComponents(1);
  bandN->SetName("b_N");
  
  vtkSmartPointer<vtkDoubleArray> bandTau = vtkSmartPointer<vtkDoubleArray>::New();
  bandTau->SetNumberOfComponents(1);
  bandTau->SetName("b_Tau");
  
  vtkSmartPointer<vtkDoubleArray> bandPress = vtkSmartPointer<vtkDoubleArray>::New();
  bandPress->SetNumberOfComponents(1);
  bandPress->SetName("b_Press");
  
  vtkSmartPointer<vtkDoubleArray> bandTensor = vtkSmartPointer<vtkDoubleArray>::New();
  bandTensor->SetNumberOfComponents(9);
  bandTensor->SetName("b_Tensor");
  
  vtkSmartPointer<vtkDoubleArray> partTensor = vtkSmartPointer<vtkDoubleArray>::New();
  partTensor->SetNumberOfComponents(9);
  partTensor->SetName("p_Tensor");
  
  vtkSmartPointer<vtkDoubleArray> bandMu = vtkSmartPointer<vtkDoubleArray>::New();
  bandMu->SetNumberOfComponents(1);
  bandMu->SetName("b_Mu");
  
  vtkSmartPointer<vtkDoubleArray> bandOmega = vtkSmartPointer<vtkDoubleArray>::New();
  bandOmega->SetNumberOfComponents(1);
  bandOmega->SetName("b_Omega");
  
  vtkSmartPointer<vtkDoubleArray> bandOmegaCoefVar = vtkSmartPointer<vtkDoubleArray>::New();
  bandOmegaCoefVar->SetNumberOfComponents(1);
  bandOmegaCoefVar->SetName("b_OmegaCoefVar");
  
  vtkSmartPointer<vtkIntArray> bandPartNum = vtkSmartPointer<vtkIntArray>::New();
  bandPartNum->SetNumberOfComponents(1);
  bandPartNum->SetName("b_PartNum");
  
  vtkSmartPointer<vtkDoubleArray> bandPartNumAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandPartNumAVG->SetNumberOfComponents(1);
  bandPartNumAVG->SetName("b_PartNumAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandVol = vtkSmartPointer<vtkDoubleArray>::New();
  bandVol->SetNumberOfComponents(1);
  bandVol->SetName("b_Vol");

  vtkSmartPointer<vtkDoubleArray> bandVolFraction = vtkSmartPointer<vtkDoubleArray>::New();
  bandVolFraction->SetNumberOfComponents(1);
  bandVolFraction->SetName("b_VolFraction");

  vtkSmartPointer<vtkDoubleArray> bandScherRate = vtkSmartPointer<vtkDoubleArray>::New();
  bandScherRate->SetNumberOfComponents(1);
  bandScherRate->SetName("b_ScherRate");

  vtkSmartPointer<vtkDoubleArray> bandContactNumAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandContactNumAVG->SetNumberOfComponents(1);
  bandContactNumAVG->SetName("b_ContactNumAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandVelLin = vtkSmartPointer<vtkDoubleArray>::New();
  bandVelLin->SetNumberOfComponents(3);
  bandVelLin->SetName("b_VelLin_rzf");
  
  vtkSmartPointer<vtkDoubleArray> bandWetContactsAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandWetContactsAVG->SetNumberOfComponents(1);
  bandWetContactsAVG->SetName("b_WetContactsAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandWetContactDistanceAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandWetContactDistanceAVG->SetNumberOfComponents(1);
  bandWetContactDistanceAVG->SetName("b_WetContactDistanceAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandVolWaterAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandVolWaterAVG->SetNumberOfComponents(1);
  bandVolWaterAVG->SetName("b_VolWaterAVG");
  
  vtkSmartPointer<vtkDoubleArray> bandVolWaterSUM = vtkSmartPointer<vtkDoubleArray>::New();
  bandVolWaterSUM->SetNumberOfComponents(1);
  bandVolWaterSUM->SetName("b_VolWaterSUM");
  
  vtkSmartPointer<vtkDoubleArray> bandDOmegaDR = vtkSmartPointer<vtkDoubleArray>::New();
  bandDOmegaDR->SetNumberOfComponents(1);
  bandDOmegaDR->SetName("b_DOmegaDR");
  
  vtkSmartPointer<vtkDoubleArray> bandOmegaNorm = vtkSmartPointer<vtkDoubleArray>::New();
  bandOmegaNorm->SetNumberOfComponents(1);
  bandOmegaNorm->SetName("b_OmegaNorm");
  
  vtkSmartPointer<vtkDoubleArray> bandGamma = vtkSmartPointer<vtkDoubleArray>::New();
  bandGamma->SetNumberOfComponents(1);
  bandGamma->SetName("b_Gamma");
  
  vtkSmartPointer<vtkDoubleArray> bandEta = vtkSmartPointer<vtkDoubleArray>::New();
  bandEta->SetNumberOfComponents(1);
  bandEta->SetName("b_Eta");
  
  vtkSmartPointer<vtkDoubleArray> bandD50 = vtkSmartPointer<vtkDoubleArray>::New();
  bandD50->SetNumberOfComponents(1);
  bandD50->SetName("b_D50");
  
  vtkSmartPointer<vtkDoubleArray> bandRadAVG = vtkSmartPointer<vtkDoubleArray>::New();
  bandRadAVG->SetNumberOfComponents(1);
  bandRadAVG->SetName("b_RadAVG");
  
  
  #ifdef ALGLIB
    vtkSmartPointer<vtkIntArray> bandShearBand = vtkSmartPointer<vtkIntArray>::New();
    bandShearBand->SetNumberOfComponents(1);
    bandShearBand->SetName("b_ShearBand");
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
    
    vtkSmartPointer<vtkIntArray> snapshot = vtkSmartPointer<vtkIntArray>::New();
    snapshot->SetNumberOfComponents(1);
    snapshot->SetName("snapshot");
    
    vtkSmartPointer<vtkFloatArray> highstressed = vtkSmartPointer<vtkFloatArray>::New();
    highstressed->SetNumberOfComponents(1);
    highstressed->SetName("highstressed");
    
    
    #ifdef ALGLIB
    vtkSmartPointer<vtkIntArray> shearband = vtkSmartPointer<vtkIntArray>::New();
    shearband->SetNumberOfComponents(1);
    shearband->SetName("shearband");
    #endif
    
    unsigned long long partId = 0;
    
    for(auto fR : _forceAll) {
      for(unsigned long long b=0; b<fR->arraySize(); b++) {
        const std::shared_ptr<particle> p1 = fR->getF(b)->part1();
        const std::shared_ptr<particle> p2 = fR->getF(b)->part2();
        
        if (not(p1->disabled()) and not(p2->disabled()) ) {
          
          partPos->InsertNextPoint(fR->getF(b)->pos1()(0), fR->getF(b)->pos1()(1), fR->getF(b)->pos1()(2));
          partPos->InsertNextPoint(fR->getF(b)->pos2()(0), fR->getF(b)->pos2()(1), fR->getF(b)->pos2()(2));
          
          vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
          line->GetPointIds()->SetId(0,(partId));
          line->GetPointIds()->SetId(1,(partId+1));
          partId+=2;
          forceCells->InsertNextCell(line);
          force->InsertNextValue(fR->getF(b)->val().norm());
          
          if (fR->getF(b)->volWater()>0) {
            wet->InsertNextValue(1);
          } else {
            wet->InsertNextValue(0);
          }
          
          if (p1->highStressedContact(p2) and p1->highStress() and p2->highStress()) {
            highstressed->InsertNextValue(fR->getF(b)->val().norm());
          } else {
            highstressed->InsertNextValue(0);
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
          snapshot->InsertNextValue(fR->getF(b)->part1()->snapshot());
        }
      }
    }
    
    vtkSmartPointer<vtkPolyData> fPd = vtkSmartPointer<vtkPolyData>::New();
    
    fPd->SetPoints(partPos);
    fPd->SetLines(forceCells);
    fPd->GetCellData()->AddArray(force);
    fPd->GetCellData()->AddArray(wet);
    fPd->GetCellData()->AddArray(snapshot);
    fPd->GetCellData()->AddArray(highstressed);
    #ifdef ALGLIB
    fPd->GetCellData()->AddArray(shearband);
    #endif
    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    
    writer->SetDataModeToAscii();
    
    #ifdef VTK6
      writer->SetInputData(fPd);
    #else
      writer->SetInput(fPd);
    #endif
    
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
            sphereSnapshot->InsertNextValue(partTemp->snapshot());
            sphereHighStress->InsertNextValue(partTemp->highStress());
            
            
            stress1->InsertNextValue(partTemp->stressSigma1());
            stress3->InsertNextValue(partTemp->stressSigma3());
            
            contactNum->InsertNextValue(partTemp->contacts());
            wetContactNum->InsertNextValue(partTemp->wetContacts());
            volWater->InsertNextValue(partTemp->volwater());
            
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
            const Eigen::Vector3d posP = Eigen::Vector3d (bandTMP->midLinedR() * cos(bandTMP->midLinedF()), bandTMP->midLinedR() * sin(bandTMP->midLinedF()), bandTMP->midLinedZ());
            pid[0] = spheresPos->InsertNextPoint(posP[0], posP[1], posP[2]);
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
          bandOmegaCoefVar->InsertNextValue(bandTMP->omegaCoefVar());
          bandPartNum->InsertNextValue(bandTMP->partNumb());
          bandPartNumAVG->InsertNextValue(bandTMP->partNumb()/bandTMP->vol());
          bandVol->InsertNextValue(bandTMP->vol());
          bandVolFraction->InsertNextValue(bandTMP->volFraction());
          bandContactNumAVG->InsertNextValue(bandTMP->contactNumAVG());
          bandScherRate->InsertNextValue(bandTMP->scherRate());
          bandMu->InsertNextValue(bandTMP->muAVG());
          bandWetContactDistanceAVG->InsertNextValue(bandTMP->wetContactDistanceAVG());
          bandWetContactsAVG->InsertNextValue(bandTMP->wetContactsAVG());
          bandVolWaterAVG->InsertNextValue(bandTMP->volwaterAVG());
          bandVolWaterSUM->InsertNextValue(bandTMP->volwaterSUM());
          bandDOmegaDR->InsertNextValue(bandTMP->dOmegadR());
          bandOmegaNorm->InsertNextValue(bandTMP->omegaNorm());
          bandGamma->InsertNextValue(bandTMP->gamma());
          bandEta->InsertNextValue(bandTMP->eta());
          bandD50->InsertNextValue(bandTMP->d50M());
          bandRadAVG->InsertNextValue(bandTMP->radAvg());
          
          double VelLin[3] = {bandTMP->vZyl()[0], bandTMP->vZyl()[1], bandTMP->vZyl()[2]};
          bandVelLin->InsertNextTupleValue(VelLin);
          #ifdef ALGLIB
            if (partTemp->shearBand()) {
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
        spheresUg->GetPointData()->AddArray(sphereSnapshot);
        spheresUg->GetPointData()->AddArray(sphereHighStress);
        spheresUg->GetPointData()->AddArray(stress1);
        spheresUg->GetPointData()->AddArray(stress3);
        spheresUg->GetPointData()->AddArray(contactNum);
        spheresUg->GetPointData()->AddArray(wetContactNum);
        spheresUg->GetPointData()->AddArray(volWater);
      }
      
      spheresUg->GetPointData()->AddArray(spheresType);
      spheresUg->GetPointData()->AddArray(bandR);
      spheresUg->GetPointData()->AddArray(bandZ);
      spheresUg->GetPointData()->AddArray(bandF);
      spheresUg->GetPointData()->AddArray(bandN);
      spheresUg->GetPointData()->AddArray(bandTau);
      spheresUg->GetPointData()->AddArray(bandPress);
      spheresUg->GetPointData()->AddArray(bandOmega);
      spheresUg->GetPointData()->AddArray(bandOmegaCoefVar);
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
      spheresUg->GetPointData()->AddArray(bandVolWaterAVG);
      spheresUg->GetPointData()->AddArray(bandVolWaterSUM);
      spheresUg->GetPointData()->AddArray(bandDOmegaDR);
      spheresUg->GetPointData()->AddArray(bandOmegaNorm);
      spheresUg->GetPointData()->AddArray(bandGamma);
      spheresUg->GetPointData()->AddArray(bandEta);
      spheresUg->GetPointData()->AddArray(bandD50);
      spheresUg->GetPointData()->AddArray(bandRadAVG);
      #ifdef ALGLIB
      spheresUg->GetPointData()->AddArray(bandShearBand);
      #endif
    }
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetDataModeToAscii();
    
    #ifdef VTK6
      writer->SetInputData(spheresUg);
    #else
      writer->SetInput(spheresUg);
    #endif
    
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
  myfileG << "042_strTensCapFR\t043_strTensCapFZ\t044_strTensCapFF\t";
  myfileG << "045_volWaterAVG\t046_volWaterSUM\t";
  myfileG << "047_dOmegadR\t048_Omega0\t049_dOmegadR\t050_Gamma\t051_OmegaCoefVar\t";
  myfileG << "052_time\t053_timeStep\t";
  myfileG << "\n";
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
    myfileG << bT->volwaterAVG()<< "\t";      // 045_volWaterAVG
    myfileG << bT->volwaterSUM()<< "\t";      // 046_volWaterSUM
    myfileG << bT->dOmegadR()<< "\t";         // 047_dOmegadR
    myfileG << bT->omega0()<< "\t";           // 048_Omega0
    myfileG << bT->dOmegadR()<< "\t";         // 049_dOmegadR
    myfileG << bT->gamma()<< "\t";            // 050_Gamma
    myfileG << bT->omegaCoefVar()<< "\t";     // 051_OmegaCoefVar
    myfileG << _cfg->_timeCur<< "\t";         // 052_time
    myfileG << _cfg->_timeStepCur<< "\t";     // 053_timeStep
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
  myfileG << "#001_id\t002_minDelta\t003_maxDelta\t004_ContNumber\t005_ContNumberAVG\t006_ForceAVG";
  
  #ifdef ALGLIB
  myfileG << "\t007_ContNumberBand\t008_ContNumberBandAVG\t009_ForceBandAVG";
  myfileG << "\t010_ContNumberBandOut\t011_ContNumberBandOutAVG\t012_ForceBandOutAVG";
  myfileG << "\t013_BandVolume\t014_OutBandVolume";
  #endif
  
  myfileG << "\t015_ContactsTotal";
  
  myfileG << "\n";
  
  std::vector <std::shared_ptr<force> >  deltas;
  std::vector <long long int>            deltasBin(bins);
  std::vector <double>                   forcesBin(bins);

  #ifdef ALGLIB
  std::vector <long long int>            deltasBinBand(bins);
  std::vector <double>                   forcesBinBand(bins);
  std::vector <long long int>            deltasBinOutBand(bins);
  std::vector <double>                   forcesBinOutBand(bins);
  const double bandVolume = _bandRow->shearBandVolume();
  const double outBandVolume = _bandRow->totalVolume() - _bandRow->shearBandVolume();
  #endif
  
  // Prepare bins for contact analyze (histogramm, force values)
  for(int x = 0; x < bins; ++x) {
    deltasBin[x] = 0;
    forcesBin[x] = 0;
  #ifdef ALGLIB
    deltasBinBand[x] = 0;
    forcesBinBand[x] = 0;
    deltasBinOutBand[x] = 0;
    forcesBinOutBand[x] = 0;
  #endif
  }
  
  double minDelta = 0.0;
  double maxDelta = 0.0;
  
  // Contact way analyze
  /////////////////////////////////////////////////
  // map for contact way analyze
  typedef boost::unordered_map<std::pair<long long, long long>, 
                       std::vector<std::pair<std::shared_ptr<snapshot>, std::shared_ptr<force>>> >
                       ContactWayMapType;
  
  ContactWayMapType contactsWay;
  long unsigned int maxVectorSize=0;
  /////////////////////////////////////////////////
  
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<force> > forces = snapshotCur->forces();
    
    for(auto f : forces) {
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
        
        // Contact way analyze
        /////////////////////////////////////////////////
        const auto contactPair = std::make_pair(std::min(f->pid1(), f->pid2()), std::max(f->pid1(), f->pid2()));
        const auto contactSnapshotForce = std::make_pair(snapshotCur, f);
        ContactWayMapType::iterator gotPair = contactsWay.find(contactPair);
        if ( gotPair == contactsWay.end() ) {
          std::vector<std::pair<std::shared_ptr<snapshot>, std::shared_ptr<force> > > newContactVector;
          newContactVector.push_back(contactSnapshotForce);
          contactsWay.insert(std::make_pair(contactPair, newContactVector));
          if (not(maxVectorSize)) maxVectorSize = 1;
        } else {
          auto tempVector = & gotPair->second;
          (*tempVector).push_back(contactSnapshotForce);
          maxVectorSize = std::max(maxVectorSize,(*tempVector).size());
        }
        /////////////////////////////////////////////////
      }
    }
  }
  
  // Contact way analyze
  /////////////////////////////////////////////////
  
  // create histogram of contact time
  std::vector<unsigned short> histContacts;
  for (long unsigned int i=0; i<maxVectorSize; i++) {histContacts.push_back(0);};
  
  for (const auto i : contactsWay ) {
    const auto tmpVector = & i.second;
    
    // Check, whether the contact was active during whole analyze-period
    unsigned int unbreakable_time = 0;
    
    if ((*tmpVector).size()> 1) {
      for (unsigned int d = 1 ; d < (*tmpVector).size(); d++) {
        const auto snap1 = (*tmpVector)[d].first;
        const auto snap0 = (*tmpVector)[d-1].first;
        if ((snap1->id() - snap0->id()) > 1) {
          histContacts[unbreakable_time]++;
          unbreakable_time = 0;
        } else {
          unbreakable_time++;
        }
      }
    }
    
    histContacts[unbreakable_time]++;
  }
  
  for (unsigned short i = histContacts.size()-1; i > 0 ; i--) {
    if (histContacts[i] == 0) {
      histContacts.erase(histContacts.begin()+i);
    } else {
      break;
    }
  }
  
  unsigned long int totContNumb = 0;
  for (const auto i : histContacts) {
    totContNumb+=i;
  }
  
  
  std::string _fileNameA;
  _fileNameA  =  _cfg->FOutput();
  _fileNameA  +=  "/contactsDuration";
  ofstream myfileA (_fileNameA.c_str());
  myfileA << "#001_stepTime\t002_ContNumber\t003_ContNumberAVG\t004_ContNumberNorm\t005_ContNumberTotal\n";
  unsigned int d = 1;
  for (const auto i : histContacts) {
    myfileA << d << "\t" << i  << "\t" << i/snapshots->size() << "\t" << double(i)/double(totContNumb) <<"\t" << totContNumb <<"\n";
    d++;
  }
  myfileA.close();
  
  /////////////////////////////////////////////////
  
  long long int allContactsTotal = 0;
  double DDelta = (maxDelta - minDelta)/bins;
  for(auto d : deltas) {
    const int binTmp = int(floor((d->deltaN()-minDelta)/DDelta));
    deltasBin[binTmp] += 1;
    forcesBin[binTmp] += d->val().norm();
    allContactsTotal++;
     #ifdef ALGLIB
      if (d->shearBand()>=0) {
        deltasBinBand[binTmp] += 1;
        forcesBinBand[binTmp] += d->val().norm();
      } else {
        deltasBinOutBand[binTmp] += 1;
        forcesBinOutBand[binTmp] += d->val().norm();
      }
    #endif
  }
  
  for(unsigned int x = 0; x < deltasBin.size(); ++x) {
    double forceTmp = 0.0;
    double forceBandTmp = 0.0;
    double forceBandOutTmp = 0.0;
    if (deltasBin[x]!=0) {
      forceTmp = forcesBin[x]/deltasBin[x];
    }
    #ifdef ALGLIB
    if (deltasBinBand[x]!=0) {
      forceBandTmp = forcesBinBand[x]/deltasBinBand[x];
    }
    if (deltasBinOutBand[x]!=0) {
      forceBandOutTmp = forcesBinOutBand[x]/deltasBinOutBand[x];
    }
    #endif
    myfileG << x << " " <<minDelta + DDelta*x  << " " <<minDelta + DDelta*(x+1) 
            << " " <<  deltasBin[x] << " " <<  deltasBin[x]/snapshots->size() 
            << " "<< forceTmp;
    #ifdef ALGLIB
    myfileG << " " <<  deltasBinBand[x] << " " <<  deltasBinBand[x]/snapshots->size() << " "<< forceBandTmp;
    myfileG << " " <<  deltasBinOutBand[x] << " " <<  deltasBinOutBand[x]/snapshots->size() << " "<< forceBandOutTmp;
    myfileG << " " <<  bandVolume << " " <<  outBandVolume;
    #endif
    myfileG << " "<< allContactsTotal;
    myfileG << "\n";
  }
  
  myfileG.close();
};

void exportclass::gnuplotContactNumberAnalyze() {  
  // Calculates, the distribution of contacts per particle
  std::string _fileNameG;
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/contactsNum";
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_ContNum\t002_PartNum\t003_PartNumNorm\t004_PartTotal";
  
  myfileG << "\n";
  
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  std::vector <long long int> contactsCalc;
  long long int partNumTotal = 0;
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    auto snapshotCur = snapshots->getSnapshot(i);
    auto particlesV = snapshotCur->particles();
    
    for (auto p : particlesV) {
      if ((p->contacts()+1) > contactsCalc.size()) {
        contactsCalc.resize(p->contacts()+1, 0);
      }
      contactsCalc[p->contacts()]++;
      partNumTotal++;
    }
  }
  
  for (unsigned int i=0; i < contactsCalc.size(); i++) {
    const double normC = (double) contactsCalc[i] / (double) partNumTotal;
    myfileG << i  << "\t" <<  contactsCalc[i] << "\t" <<  normC
                  << "\t" <<  partNumTotal << "\n";
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
    
    for(auto f : forces) {
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
    
    for(auto f : forces) {
      std::shared_ptr<contactFollow> tmpCntFolw = std::make_shared<contactFollow>(f, snapshotCur);
      cntFolw.push_back(tmpCntFolw);
    }
  }
  
  std::sort (cntFolw.begin(), cntFolw.end(), sortContactFollow);
  
  for(auto tmpCntFolw : cntFolw) {
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
    
    for(auto p : particles) {
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
    
    
    for(auto f : forces) {
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
      // const double ThetaId = ThetaI*dAngle + d2Ang;    // Not used for the moment
      
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
  
  myfile2 << "#001_Id\t002_time\t003_torque\t004_kinEnergy\t005_potEnergy\n";
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    myfile2                                                             << i << "\t";  // 001_snapId
    myfile2                            << snapshotCur->timeStep()*_cfg->dT() << "\t";  // 002_time
    myfile2 << snapshotCur->torque(_cfg->get_o(), _cfg->get_c(), _cfg->tR()) << "\t";  // 003_torque
    myfile2 << snapshotCur->kinEnergy(_cfg->tC()) << "\t";                             // 004_kinEnergy
    myfile2 << snapshotCur->potEnergy(_cfg->tC()) << "\n";                             // 005_potEnergy
  }
  myfile2.close();
}

void exportclass::forceChain() {
  if (not(_noForces)) {
    std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
    for(unsigned int i=0; i<snapshots->size(); i++) {
      std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
      snapshotCur->forceChainRet();
    }
  }
}


void exportclass::gnuplotWetParticles(int bins) {  
  // Calculates, how much particles have liquid, sorted by water volume
  std::string _fileNameG;
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/wetParticles";
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#001_id\t002_minWat\t003_maxWat\t004_partNumber\t005_relPartNumber";
  
  myfileG << "\n";
  
  std::vector <long long int>            deltasBin(bins);
  
  for(int x = 0; x < bins; ++x) {
    deltasBin[x] = 0;
  }
  
  double minWat = 0.0;
  double maxWat = 0.0;
    
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<particle> > particles = snapshotCur->particles();
    
    for(auto p : particles) {
      if (p->disabled() or p->volwater()<=0) {continue;}
      
      if  (minWat == maxWat and minWat == 0 ) {
        minWat = p->relativeWaterVol();
        maxWat = p->relativeWaterVol();
      }
      minWat = std::min(minWat, p->relativeWaterVol());
      maxWat = std::max(maxWat, p->relativeWaterVol());
    }
  }
  
  double DDelta = (maxWat - minWat)/bins;
  long long int partNumb = 0;
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    std::vector <std::shared_ptr<particle> > particles = snapshotCur->particles();
    
    for(auto p : particles) {
      if (p->disabled() or p->volwater()<=0) {continue;}
      
      deltasBin[int(floor((p->relativeWaterVol()-minWat)/DDelta))] += 1;
      partNumb++;
    }
  }
  
  for(unsigned int x = 0; x < deltasBin.size(); ++x) {
    myfileG << x << "\t" << minWat + DDelta*x << "\t" <<  minWat + DDelta*(x+1) 
            <<      "\t" << deltasBin[x]      << "\t" <<  deltasBin[x]/double(partNumb) ;
    myfileG << "\n";
  }
  
  myfileG.close();
};
