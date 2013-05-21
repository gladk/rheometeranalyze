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
    spheresUg->GetPointData()->AddArray(bandVelLin);
    spheresUg->GetPointData()->AddArray(bandTensor);
    spheresUg->GetPointData()->AddArray(partTensor);
    spheresUg->GetPointData()->AddArray(posZyl);
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
  myfileG << "031_ShearBand\t032_WetContactsAVG\t033_WetContactsDistAVG\t \n";
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
    myfileG << " \n"; 
  }
        
  myfileG.close();
};

void exportclass::Utwente()  {
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  //unsigned int leadingZerosLen = static_cast <unsigned int> (log10 (snapshots->size()) + 1);
  unsigned int leadingZerosLen = 4;
  
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
              << "0.0" << " "                               //Theta, should be fixed.
              << f->valN() << " " <<  "0.0" << " "          //Normal and tangential components of the force
              << f->nVec()(0) << " " << f->nVec()(1) << " " << f->nVec()(2) << " " 
              << "0.0" << " " << "0.0" << " " << "0.0" << " " 
              << std::endl;
      } else {
        Fstat << timeTmp << " " << pid2T << " " << pid1T << " "
              << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
              << f->deltaN() << " "
              << "0.0" << " "                               //Theta, should be fixed.
              << f->valN() << " " <<  "0.0" << " "          //Normal and tangential components of the force
              << -f->nVec()(0) << " " << -f->nVec()(1) << " " << -f->nVec()(2) << " " 
              << "0.0" << " " << "0.0" << " " << "0.0" << " " 
              << std::endl;
      }
      if ((pid1T>=0) and (pid2T>=0)) {
        Fstat << timeTmp << " " << pid2T << " " << pid1T << " "
              << f->cP()(0) << " " << f->cP()(1) << " " << f->cP()(2) << " " 
              << f->deltaN() << " "
              << "0.0" << " "                               //Theta, should be fixed.
              << f->valN() << " " <<  "0.0" << " "          //Normal and tangential components of the force
              << -f->nVec()(0) << " " << -f->nVec()(1) << " " << -f->nVec()(2) << " " 
              << "0.0" << " " << "0.0" << " " << "0.0" << " " 
              << std::endl;
      }
    }
    
    Fstat.close();
  }
}
