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

#include "rheometer.h"
#include "snapshot.h"

rheometer::rheometer(std::shared_ptr<configopt> cfg) {
  _cfg = cfg;
  _particleNum = 0;
  _forceNum = 0;
  
  loadParticles();
  
  //Create bands
  std::shared_ptr <bandRow> bandRowTMP (new bandRow(_cfg, _particleAll,  _forceRow));
  std::shared_ptr <bandRow> _bandRow = bandRowTMP;
  
  
  std::shared_ptr <exportclass> exp (new exportclass(_cfg, _bandRow));
  
  if (_cfg->Vtk()) exp->VTK();
  if (_cfg->Utwente()) exp->Utwente();
  
  exp->gnuplotSchearRate();
};

void rheometer::loadParticles() {
  unsigned int partNumbCounter  = 1;
  std::shared_ptr<snapshotRow> snapshots = _cfg->snapshot();
  
  for(unsigned int i=0; i<snapshots->size(); i++) {
    
    std::shared_ptr<snapshot> snapshotCur = snapshots->getSnapshot(i);
    
    std::ifstream _file;
    _file.open(snapshotCur->getParticleFile().string());
    
    std::string   line;
    int curLine = 1;
    unsigned long long maxId = 0;
    std::vector <std::shared_ptr<particle> > tmpPartVector;
  
    while(std::getline(_file, line)) {
      std::stringstream linestream(line);
      std::string data;
      
      int valInt;
      double valD;
      double pR, pM, pD;
      int pT;
      unsigned long long pId;
      Eigen::Vector3f pC, pV, pO;
      if (curLine>=_cfg->nDat()) {
        for (int i=1; i<=_cfg->maxC(); i++) {
          if (i==_cfg->cId()) {
            linestream >> pId;
          } else if (i==_cfg->cT()) {
            linestream >> pT;
            //std::cerr<<pT<<std::endl;
          } else if (i==_cfg->cC()) {
            linestream >> pC[0];
            linestream >> pC[1];
            linestream >> pC[2];
            i+=2;
            //std::cerr<<pC<<std::endl<<std::endl;
          } else if (i==_cfg->cV()) {
            linestream >> pV[0];
            linestream >> pV[1];
            linestream >> pV[2];
            i+=2;
            //std::cerr<<pV<<std::endl<<std::endl;
          } else if (i==_cfg->cO()) {
            linestream >> pO[0];
            linestream >> pO[1];
            linestream >> pO[2];
            i+=2;
            //std::cerr<<pO<<std::endl<<std::endl;
          } else if (i==_cfg->cR()) {
            linestream >> pR;
            //std::cerr<<pR<<std::endl;
          } else if (i==_cfg->cM()) {
            linestream >> pM;
            //std::cerr<<pM<<std::endl;
          } else if (i==_cfg->cD()) {
            linestream >> pD;
            //std::cerr<<pD<<std::endl;
          } else {
            linestream >> valD;
          }
        }
        
        if (((_cfg->tC()>=0) and (pT == _cfg->tC())) or (_cfg->tC()<0)) {
          maxId = max(pId, maxId);
          std::shared_ptr<particle> tmpParticle ( new particle (pId, pT, partNumbCounter-1, pR, pM, pD, pC,pV, pO));
          tmpPartVector.push_back(tmpParticle);
          snapshotCur->addParticle(tmpParticle);
        }
  
      } else if (curLine == _cfg->nAt()) {
        linestream >> valInt;
        std::cerr<<"File "<<partNumbCounter<<"/"<<snapshots->size()<<std::endl<<"Expected particles "<<valInt;
      } else if (curLine == _cfg->nPSt()) {
        linestream >> valInt;
        snapshotCur->setTimeStep(valInt);
      }
      curLine++;
    };
    std::shared_ptr<particleRow> particleTMP ( new particleRow(maxId+1));
    _particleAll.push_back(particleTMP);
    
    unsigned int partNumbTMP  = _particleAll.size()-1;
    BOOST_FOREACH( std::shared_ptr<particle> p, tmpPartVector) {
       _particleAll[partNumbTMP]->addP(p);
    }
    std::cerr<<"; "<<_particleAll[partNumbTMP]->elementsNum()<<" particles added."<<std::endl;
    _particleNum+=_particleAll[partNumbTMP]->elementsNum();
    partNumbCounter++;
    
    this->loadForces(snapshotCur);
  
  }
  snapshots->sortRow();
  std::cerr<<"The total number of added particles is "<<_particleNum<<std::endl;
  std::cerr<<"The total number of added forces is "<<_forceNum<<std::endl;
};


void rheometer::loadForces(std::shared_ptr<snapshot> loadSnap) {
    std::ifstream _file;
    _file.open(loadSnap->getForceFile().string());
    
    std::string   line;
    int curLine = 1;
    
    //std::vector <std::shared_ptr<particle> > tmpPartVector;
    
    int valInt;
    double valD;
    
    std::shared_ptr<forceRow> forceTMP ( new forceRow());
    _forceRow.push_back(forceTMP);
    unsigned int forceRowNumbTMP  = _forceRow.size()-1;
    
    
    while(std::getline(_file, line)) {
      std::stringstream linestream(line);
      std::string data;
      
      unsigned long long pid1, pid2;
      Eigen::Vector3f pos1, pos2, val;
    
      if (curLine>=_cfg->fDat()) {
        for (int i=1; i<=_cfg->maxCF(); i++) {
          if (i==_cfg->cPos1ID()) {
            linestream >> pid1;
          } else if (i==_cfg->cPos2ID()) {
            linestream >> pid2;
          } else if (i==_cfg->cPos1()) {
            linestream >> pos1[0];
            linestream >> pos1[1];
            linestream >> pos1[2];
            i+=2;
            //std::cerr<<"pos1 " << pos1<<std::endl<<std::endl;
          } else if (i==_cfg->cPos2()) {
            linestream >> pos2[0];
            linestream >> pos2[1];
            linestream >> pos2[2];
            i+=2;
            //std::cerr<<"pos2 " << pos2<<std::endl<<std::endl;
          } else if (i==_cfg->cForc()) {
            linestream >> val[0];
            linestream >> val[1];
            linestream >> val[2];
            i+=2;
            //std::cerr<<"val " << val<<std::endl<<std::endl;
          } else {
            linestream >> valD;
          }
        };
        
        if (_particleAll[forceRowNumbTMP]->particleReal(pid1) and _particleAll[forceRowNumbTMP]->particleReal(pid2)){
          std::shared_ptr<force> tmpForce ( new force (_particleAll[forceRowNumbTMP]->getP(pid1), _particleAll[forceRowNumbTMP]->getP(pid2), forceRowNumbTMP, pos1, pos2, val));
          _forceRow[forceRowNumbTMP]->addF(tmpForce);
          loadSnap->addForce(tmpForce);
        }
  
      } else if (curLine == _cfg->fAt()) {
        linestream >> valInt;
         std::cerr<<"Expected forces "<<valInt;
      } else if (curLine == _cfg->nFSt()) {
        linestream >> valInt;
        unsigned int valIntUn = valInt;
        if (valIntUn != loadSnap->timeStep()) {
          std::cerr << "Timestep of force and particle files is not the same!\n"; 
          exit (EXIT_FAILURE);
        }
      }
      curLine++;
    };
    std::cerr<<"; "<<_forceRow[forceRowNumbTMP]->elementsNum()<<" forces added."<<std::endl<<std::endl;
    _forceNum+=_forceRow[forceRowNumbTMP]->elementsNum();
};
