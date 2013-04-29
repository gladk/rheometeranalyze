/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "export.h"
#include <boost/foreach.hpp>

exportclass::exportclass(std::shared_ptr<configopt> cfg, std::shared_ptr <bandRow> bandAll) {
  _cfg = cfg;
  _bandRow = bandAll;
};

void exportclass::VTK() {
  ofstream fileVTK;
  std::string _fileNameVTU;

  _fileNameVTU =  _cfg->FOutput();
  _fileNameVTU += "/output.vtu";

  
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
  
  vtkSmartPointer<vtkDoubleArray> bandVelLinDf = vtkSmartPointer<vtkDoubleArray>::New();
  bandVelLinDf->SetNumberOfComponents(1);
  bandVelLinDf->SetName("bandVelLinDf");
  
  #ifdef ALGLIB
    vtkSmartPointer<vtkIntArray> bandShearBand = vtkSmartPointer<vtkIntArray>::New();
    bandShearBand->SetNumberOfComponents(1);
    bandShearBand->SetName("bandShearBand");
  #endif
  
  vtkSmartPointer<vtkUnstructuredGrid> spheresUg = vtkSmartPointer<vtkUnstructuredGrid>::New();
  
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bandTMP = _bandRow->getBand(b);
    
    for (int z = 0; z<bandTMP->partNumb(); z++) {
      std::shared_ptr<particle> partTemp = bandTMP->getPart(z);
      if (not(partTemp->disabled())) {
        vtkIdType pid[1];
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
        
        bandR->InsertNextValue(bandTMP->idR());
        bandZ->InsertNextValue(bandTMP->idZ());
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
        bandVelLinDf->InsertNextValue(bandTMP->vDf());
        #ifdef ALGLIB
          if (bandTMP->shearBand()) {
            bandShearBand->InsertNextValue(1);
          } else {
            bandShearBand->InsertNextValue(0);
          }
        #endif
        
        spheresCells->InsertNextCell(1,pid);
      }
    }
    
    
    spheresUg->SetPoints(spheresPos);
    spheresUg->SetCells(VTK_VERTEX, spheresCells);
    spheresUg->GetPointData()->AddArray(radii);
    spheresUg->GetPointData()->AddArray(mass);
    spheresUg->GetPointData()->AddArray(density);
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
    spheresUg->GetPointData()->AddArray(bandOmega);
    spheresUg->GetPointData()->AddArray(bandPartNum);
    spheresUg->GetPointData()->AddArray(bandPartNumAVG);
    spheresUg->GetPointData()->AddArray(bandVol);
    spheresUg->GetPointData()->AddArray(bandVolFraction);
    spheresUg->GetPointData()->AddArray(bandContactNumAVG);
    spheresUg->GetPointData()->AddArray(bandScherRate);
    spheresUg->GetPointData()->AddArray(bandMu);
    spheresUg->GetPointData()->AddArray(bandVelLinDf);
    spheresUg->GetPointData()->AddArray(bandTensor);
    spheresUg->GetPointData()->AddArray(partTensor);
    spheresUg->GetPointData()->AddArray(posZyl);
    #ifdef ALGLIB
    spheresUg->GetPointData()->AddArray(bandShearBand);
    #endif
  }
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetInput(spheresUg);
  
  writer->SetFileName(_fileNameVTU.c_str());
  writer->Write();
};

void exportclass::gnuplotSchearRate() {  
  
  std::string _fileNameG1, _fileNameG2, _fileNameG3, _fileNameG4, _fileNameG5, _fileNameG;
  _fileNameG1 =  _cfg->FOutput();
  _fileNameG1 += "/gnuplot_shearrate.txt";
  
  _fileNameG2 =  _cfg->FOutput();
  _fileNameG2 += "/gnuplot_scherrate2.txt";
  
  _fileNameG3 =  _cfg->FOutput();
  _fileNameG3 += "/gnuplot_vel1_c.txt";
  
  _fileNameG4 =  _cfg->FOutput();
  _fileNameG4 += "/shearband";
  
  _fileNameG5 =  _cfg->FOutput();
  _fileNameG5 += "/gnuplot_vel2.txt";
  
  _fileNameG  =  _cfg->FOutput();
  _fileNameG  +=  "/gnuplot_daten";
  
  ofstream myfile1 (_fileNameG1.c_str());
  ofstream myfile11 (_fileNameG2.c_str());
  ofstream myfile4 (_fileNameG4.c_str());
  ofstream myfile5 (_fileNameG5.c_str());
  
  //HACK==================================
  std::string myfile001_name = _fileNameG2; myfile001_name += "_001"; ofstream myfile001 (myfile001_name.c_str());
  std::string myfile002_name = _fileNameG2; myfile002_name += "_002"; ofstream myfile002 (myfile002_name.c_str());
  std::string myfile004_name = _fileNameG2; myfile004_name += "_004"; ofstream myfile004 (myfile004_name.c_str());
  std::string myfile006_name = _fileNameG2; myfile006_name += "_006"; ofstream myfile006 (myfile006_name.c_str());
  std::string myfile008_name = _fileNameG2; myfile008_name += "_008"; ofstream myfile008 (myfile008_name.c_str());
  std::string myfile016_name = _fileNameG2; myfile016_name += "_016"; ofstream myfile016 (myfile016_name.c_str());
  std::string myfile100_name = _fileNameG2; myfile100_name += "_100"; ofstream myfile100 (myfile100_name.c_str());
  //HACK==================================
  
  
  if (myfile1.is_open() and myfile11.is_open()) {
    myfile1 << "P\tSigmaD\tmu" << std::endl;
    myfile11 << "PressGlob\tTauGlob\tScherrate\tScherrateI\tEta\tXPos\tZPos" << std::endl;
    for(unsigned int b=0; b<_bandRow->size(); b++) {
      std::shared_ptr<band> bandTMP = _bandRow->getBand(b);
      if (bandTMP->partNumb()>0) {
        double p = bandTMP->PressAVG();
        double mu = bandTMP->muAVG();
        
        bandTMP->midLinedZ();
        
        if (bandTMP->scherRate() > 0.0) {
          myfile1 << p << "\t"<< p << "\t"<< mu << std::endl;
          myfile11 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
        }

        //HACK==================================
        if (bandTMP->scherRate() >= 0.0) {
          int shearBandT;
          if (bandTMP->shearBand()){
            shearBandT = 1;
          } else {
            shearBandT = 0;
          }
          
          if (bandTMP->scherRate() <= 0.01) {
            myfile001 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else if (bandTMP->scherRate() <= 0.02) {
            myfile002 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else if (bandTMP->scherRate() <= 0.04) {
            myfile004 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else if (bandTMP->scherRate() <= 0.06) {
            myfile006 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else if (bandTMP->scherRate() <= 0.08) {
            myfile008 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else if (bandTMP->scherRate() <= 0.16) {
            myfile016 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          } else {
            myfile100 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << "\t"<< shearBandT << std::endl;
          }
        }
        //HACK==================================
        
      }
    }
  }
  
  #ifdef ALGLIB
  if (myfile4.is_open()) {
    myfile4 << "RZ\tW\tH" << std::endl;
    for(unsigned int b=0; b<_bandRow->shearBandSize(); b++) {
      const Eigen::Vector3d tmpVal = _bandRow->shearBand(b);
      myfile4 << tmpVal[0] << "\t"<< tmpVal[1] << "\t"<< tmpVal[2] << std::endl;
    }
  }
  #endif
  
  myfile1.close();
  myfile11.close();
  myfile4.close();

  myfile001.close();
  myfile002.close();
  myfile004.close();
  myfile006.close();
  myfile008.close();
  myfile016.close();
  myfile100.close();
  
  ofstream myfile2 (_fileNameG3.c_str());
  if (myfile2.is_open()) {
    myfile2 << "y_durch_d\ty_durch_d^2\tVabs\tV_durch_V0\tV_Deviation\tVolFraction\td\tScherrate" << std::endl;
    double V0 = 0;
    for(int R=0; R<_cfg->SecRadial(); R++) {
      double rAVG = 0;
      double midLinedR = 0;
      double y_durch_d = 0;
      double V = 0;
      double VolFract = 0;
      double VStDev = 0;
      double ScherR = 0;
      
      for(int Z=0; Z<_cfg->SecZ(); Z++) {
        std::shared_ptr<band> bandTMP = _bandRow->getBand(R, Z);
        rAVG += bandTMP->radAvg();
        midLinedR += bandTMP->midLinedR();
        V += bandTMP->omega();
        VolFract += bandTMP->volFraction();
        VStDev += bandTMP->omegaStDev();
        ScherR += bandTMP->scherRate();
      }
      rAVG /=_cfg->SecRadial();
      V /=_cfg->SecRadial();
      VolFract /=_cfg->SecRadial();
      VStDev /=_cfg->SecRadial();
      ScherR /=_cfg->SecRadial();
      if (R == 0) {V0 = V;};
      midLinedR = midLinedR/_cfg->SecRadial()  - _cfg->Din()/2.0;
      y_durch_d = midLinedR/(rAVG*2.0);
      myfile2 << y_durch_d<< "\t"<< y_durch_d*y_durch_d << "\t"<< V << "\t"<< V/V0 << "\t"<< VStDev << "\t"<< VolFract << "\t"<< rAVG*2.0<< "\t"<< ScherR << std::endl;
    }
  }  
  myfile2.close();
  
  if (myfile5.is_open()) {
    myfile5 << "#R\t";
    for(unsigned int h=0; h<_cfg->SecZ(); h++) {
      myfile5 << "H="<<(_bandRow->getBand(0,h))->midLinedZ()<<"\t";
    }
    myfile5 << std::endl;
    
    double omega0 = _bandRow->getBand(_cfg->SecRadial()-1,0)->omega(); 
    for(unsigned int r=0; r<_cfg->SecRadial(); r++) {
      myfile5 << (_bandRow->getBand(r,0))->midLinedR()<<"\t";
      for(unsigned int h=0; h<_cfg->SecZ(); h++) {
        double omega = _bandRow->getBand(r,h)->omega(); 
        double valT = omega/omega0;
        myfile5 << (valT)<<"\t";
      }
      myfile5 << std::endl;
    }
  }
  myfile5.close();
  
  ofstream myfileG (_fileNameG.c_str());
  myfileG << "#000_id\t001_r\t002_z\t003_rPos\t004_zPos\t005_H\t006_W\t007_PartN\t";
  myfileG << "008_vol\t009_volFract\t010_contactNum\t011_vZyldR\t012_vZyldZ\t013_vZyldF\t";
  myfileG << "014_ShearRate\t015_I\t016_Eta\t017_omega\t018_press\t019_tau\t020_mu\t";
  myfileG << "021_strTensRR\t022_strTensRZ\t023_strTensRF\t";
  myfileG << "024_strTensZR\t025_strTensZZ\t026_strTensZF\t";
  myfileG << "027_strTensFR\t028_strTensFZ\t029_strTensFF\t";
  myfileG << "030_ShearBand\t \n";
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bT = _bandRow->getBand(b);
    myfileG << bT->id() << "\t";           // 000_id
    myfileG << bT->idR() << "\t";          // 001_r
    myfileG << bT->idZ() << "\t";          // 002_z
    myfileG << bT->midLinedR() << "\t";    // 003_rPos
    myfileG << bT->midLinedZ() << "\t";    // 004_zPos
    myfileG << bT->dZ() << "\t";           // 005_H
    myfileG << bT->dR() << "\t";           // 006_W
    myfileG << bT->partNumb() << "\t";     // 007_PartN
    myfileG << bT->vol() << "\t";          // 008_vol
    myfileG << bT->volFraction() << "\t";  // 009_volFract
    myfileG << bT->contactNumAVG() << "\t";// 010_contactNum
    myfileG << bT->vZyl()[0]<< "\t";       // 011_vZyldR
    myfileG << bT->vZyl()[1]<< "\t";       // 012_vZyldZ
    myfileG << bT->vZyl()[2]<< "\t";       // 013_vZyldF
    myfileG << bT->scherRate()<< "\t";     // 014_ShearRate
    myfileG << bT->I()<< "\t";             // 015_I
    myfileG << bT->eta()<< "\t";           // 016_Eta
    myfileG << bT->omega()<< "\t";         // 017_omega
    myfileG << bT->press()<< "\t";         // 018_press
    myfileG << bT->tau()<< "\t";           // 019_tau
    myfileG << bT->muAVG()<< "\t";         // 020_mu
    myfileG << bT->TensorAVG()(0)<< "\t";  // 021_strTensRR
    myfileG << bT->TensorAVG()(1)<< "\t";  // 022_strTensRZ
    myfileG << bT->TensorAVG()(2)<< "\t";  // 023_strTensRF
    myfileG << bT->TensorAVG()(3)<< "\t";  // 024_strTensZR
    myfileG << bT->TensorAVG()(4)<< "\t";  // 025_strTensZZ
    myfileG << bT->TensorAVG()(5)<< "\t";  // 026_strTensZF
    myfileG << bT->TensorAVG()(6)<< "\t";  // 027_strTensFR
    myfileG << bT->TensorAVG()(7)<< "\t";  // 028_strTensFZ
    myfileG << bT->TensorAVG()(8)<< "\t";  // 029_strTensFF
    myfileG << bT->shearBand()<< "\t";     // 030_ShearBand - True (1) or False (0)
    myfileG << " \n"; 
  }
        
  myfileG.close();
};

void exportclass::Utwente()  {
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  unsigned int leadingZerosLen = static_cast <unsigned int> (log10 (snapshots->size()) + 1);
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    
    stringstream ss;
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
          << p->rad() << std::endl;
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
    Fstat<<"# TimeStep " << snapshotCur->timeStep() << std::endl;
    Fstat<<"# Time " << timeTmp << std::endl;
    Fstat<<"#  " << std::endl;
    
    
    BOOST_FOREACH(std::shared_ptr<force> f, forces) {
      Fstat << timeTmp << " " << f->pid1() << " " << f->pid2() << " "
            << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
            << f->deltaN() << " "
            << "0.0" << " "                                      //Theta, should be fixed.
            << f->valN() << " " <<  "0.0" << " "          //Normal and tangential components of the force
            << f->nVec()(0) << " " << f->nVec()(1) << " " << f->nVec()(2) << " " 
            << "0.0" << " " << "0.0" << " " << "0.0" << " " 
            << std::endl;
    }
    
    Fstat.close();
  }
}
