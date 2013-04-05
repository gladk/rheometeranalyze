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
  
  
  _fileNameVTU =  _cfg->FOutput();
  _fileNameVTU += "/output.vtu";
  
  _fileNameG1 =  _cfg->FOutput();
  _fileNameG1 += "/gnuplot_shearrate.txt";
  
  _fileNameG2 =  _cfg->FOutput();
  _fileNameG2 += "/gnuplot_scherrate2.txt";
  
  _fileNameG3 =  _cfg->FOutput();
  _fileNameG3 += "/gnuplot_vel1_c.txt";
  
  _fileNameG4 =  _cfg->FOutput();
  _fileNameG4 += "/gnuplot_shearZone.txt";
  
  
};

void exportclass::VTK() {
  
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
  
  vtkSmartPointer<vtkDoubleArray> bandGlobMu = vtkSmartPointer<vtkDoubleArray>::New();
  bandGlobMu->SetNumberOfComponents(1);
  bandGlobMu->SetName("bandGlobMu");
  
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
  
  vtkSmartPointer<vtkDoubleArray> bandVelLinDf = vtkSmartPointer<vtkDoubleArray>::New();
  bandVelLinDf->SetNumberOfComponents(1);
  bandVelLinDf->SetName("bandVelLinDf");
  
  
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
        bandGlobMu->InsertNextValue(bandTMP->muGlobAVG());
        bandVelLinDf->InsertNextValue(bandTMP->vDf());
        
        spheresCells->InsertNextCell(1,pid);
      }
    }
    
    
    spheresUg->SetPoints(spheresPos);
    spheresUg->SetCells(VTK_VERTEX, spheresCells);
    spheresUg->GetPointData()->AddArray(radii);
    spheresUg->GetPointData()->AddArray(mass);
    spheresUg->GetPointData()->AddArray(density);
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
    spheresUg->GetPointData()->AddArray(bandGlobMu);
    spheresUg->GetPointData()->AddArray(bandVelLinDf);
    
  }
  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetDataModeToAscii();
  writer->SetInput(spheresUg);
  
  writer->SetFileName(_fileNameVTU.c_str());
  writer->Write();
};

void exportclass::gnuplotSchearRate() {  
  ofstream myfile1 (_fileNameG1.c_str());
  ofstream myfile11 (_fileNameG2.c_str());
  ofstream myfile4 (_fileNameG4.c_str());
  
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
        double p = bandTMP->localPress();
        double SigmD = bandTMP->dLocalAvg();
        double mu = bandTMP->muLocalAVG();
        
        
        bandTMP->midLinedZ();
        
        if (bandTMP->scherRate() > 0.0) {
          myfile1 << p << "\t"<< SigmD << "\t"<< mu << std::endl;
          myfile11 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
        }

        //HACK==================================
        if (bandTMP->scherRate() >= 0.0) {
          if (bandTMP->scherRate() <= 0.01) {
            myfile001 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else if (bandTMP->scherRate() <= 0.02) {
            myfile002 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else if (bandTMP->scherRate() <= 0.04) {
            myfile004 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else if (bandTMP->scherRate() <= 0.06) {
            myfile006 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else if (bandTMP->scherRate() <= 0.08) {
            myfile008 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else if (bandTMP->scherRate() <= 0.16) {
            myfile016 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          } else {
            myfile100 << bandTMP->press() << "\t"<< bandTMP->tau() << "\t"<< bandTMP->scherRate() << "\t"<< bandTMP->I()  << "\t"<< bandTMP->eta() << "\t"<< bandTMP->midLinedR() <<  "\t"<< bandTMP->midLinedZ() << std::endl;
          }
        }
        //HACK==================================
        
      }
    }
  }
  
  if (myfile4.is_open()) {
    myfile4 << "RPOS\tZPOS\tW" << std::endl;
    for(unsigned int b=0; b<_bandRow->getBandShearZonesSize(); b++) {
      std::shared_ptr<bandShearZone> bandSZ = _bandRow->getBandShearZones(b);
      myfile4 << bandSZ->RPOS() << "\t"<< bandSZ->ZPOS() << "\t"<< bandSZ->W() << std::endl;
    }
  }
    
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
    C3d << particles.size() << "\t" << timeTmp
        << "\t" << -_cfg->Dout()/2.0 << "\t"<<-_cfg->Dout()/2.0 << "\t0.0\t"
        << _cfg->Dout()/2.0 << "\t" << _cfg->Dout()/2.0 << "\t" << _cfg->H()
        << std::endl;
    
    BOOST_FOREACH(std::shared_ptr<particle> p, particles) {
      C3d << p->c()(0) << "\t" << p->c()(1) << "\t"<< p->c()(2) << "\t" 
          << p->v()(0) << "\t" << p->v()(1) << "\t"<< p->v()(2) << "\t"
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
      Fstat << timeTmp << "\t" << f->pid1() << "\t" << f->pid2() << "\t"
            << f->cP()(0) << "\t" << f->cP()(1) << "\t" << f->cP()(2) << "\t" 
            << f->deltaN() << "\t"
            << "0.0" << "\t"                                      //Theta, should be fixed.
            << f->valN() << "\t" <<  "0.0" << "\t"          //Normal and tangential components of the force
            << f->nVec()(0) << "\t" << f->nVec()(1) << "\t" << f->nVec()(2) << "\t" 
            << "0.0" << "\t" << "0.0" << "\t" << "0.0" << "\t" 
            << std::endl;
    }
    
    Fstat.close();
    /*
    C3d << particles.size() << "\t" << snapshotCur->timeStep()*_cfg->dT()
        << "\t" << -_cfg->Dout()/2.0 << "\t"<<-_cfg->Dout()/2.0 << "\t0.0\t"
        << _cfg->Dout()/2.0 << "\t" << _cfg->Dout()/2.0 << "\t" << _cfg->H()
        << std::endl;
    
    BOOST_FOREACH(std::shared_ptr<particle> p, particles) {
      C3d << p->c()(0) << "\t" << p->c()(1) << "\t"<< p->c()(2) << "\t" 
          << p->v()(0) << "\t" << p->v()(1) << "\t"<< p->v()(2) << "\t"
          << p->rad() << std::endl;
    };
    
    C3d << "\t"<< std::endl;
    C3d.close();
    */ 
    //================================================
  }
}
