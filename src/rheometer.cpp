/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014 Anton Gladky <gladky.anton@gmail.com>

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
  
  std::shared_ptr <exportclass> exp (new exportclass(_cfg, _bandRow, _forceRow));
  
  // ==========================ForceChain
  
  exp->forceChain();
  
  // ==========================ForceChain
  
  
  if (_cfg->Vtk()) exp->VTK();
  if (_cfg->Utwente()) exp->Utwente();
  
  exp->gnuplotSchearRate();
  exp->torque();
  if (_cfg->contact()) {
    exp->gnuplotContactAnalyze(100);
    exp->gnuplotContactWet();
  }
  
  if (_cfg->followContact()) {
    exp->gnuplotContactFollow();
  }
  
  if (_cfg->intOri()>0) {
    exp->intOri();
  }
  
};

void rheometer::loadParticles() {
  namespace src = boost::log::sources;
  namespace logging = boost::log;
  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;
  
  
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
      double pR, pM, pD, volWater=0.0;
      int pT;
      unsigned long long pId;
      Eigen::Vector3d pC, pV, pO;
      
      if (curLine>=_cfg->nDat()) {
        for (int i=1; i<=_cfg->maxC(); i++) {
          if (i==_cfg->cId()) {
            linestream >> pId;
          } else if (i==_cfg->cT()) {
            linestream >> pT;
            //std::cerr<<pT;
          } else if (i==_cfg->cC()) {
            linestream >> pC[0];
            linestream >> pC[1];
            linestream >> pC[2];
            i+=2;
            //std::cerr<<pC<<std::endl;
          } else if (i==_cfg->cV()) {
            linestream >> pV[0];
            linestream >> pV[1];
            linestream >> pV[2];
            i+=2;
            //std::cerr<<pV<<std::endl;
          } else if (i==_cfg->cO()) {
            linestream >> pO[0];
            linestream >> pO[1];
            linestream >> pO[2];
            i+=2;
            //std::cerr<<pO<<std::endl;
          } else if (i==_cfg->cR()) {
            linestream >> pR;
            //std::cerr<<pR;
          } else if (i==_cfg->cM()) {
            linestream >> pM;
            //std::cerr<<pM;
          } else if (i==_cfg->cD()) {
            linestream >> pD;
            //std::cerr<<pD;
          } else if (i==_cfg->cVolWaterP()) {
            linestream >> volWater;
            //std::cerr<<"VolWaterP " << VolWaterP<<std::endl;
          }  else {
            linestream >> valD;
          }
        }
        
        if (((_cfg->tC()>=0) and (pT == _cfg->tC())) or (_cfg->tC()<0)) {
          maxId = max(pId, maxId);
          std::shared_ptr<particle> tmpParticle ( new particle (pId, pT, partNumbCounter-1, pR, pM, pD, pC,pV, pO, volWater));
          if (_cfg->intOri() > 0) tmpParticle->createIntOri(_cfg->intOri());
          tmpPartVector.push_back(tmpParticle);
          snapshotCur->addParticle(tmpParticle);
        }
  
      } else if (curLine == _cfg->nAt()) {
        linestream >> valInt;
        BOOST_LOG_SEV(lg, info)<<"File "<<partNumbCounter<<"/"<<snapshots->size()<< " (" << snapshotCur->getParticleFile() << "); " <<"Expected particles "<<valInt;
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
    BOOST_LOG_SEV(lg, info)<<_particleAll[partNumbTMP]->elementsNum()<<" particles added.";
    _particleNum+=_particleAll[partNumbTMP]->elementsNum();
    partNumbCounter++;
    
    this->loadForces(snapshotCur);
  
  }
  snapshots->sortRow();
  BOOST_LOG_SEV(lg, info)<<"The total number of added particles is "<<_particleNum;
  BOOST_LOG_SEV(lg, info)<<"The total number of added forces is "<<_forceNum;
};


void rheometer::loadForces(std::shared_ptr<snapshot> loadSnap) {
  namespace src = boost::log::sources;
  namespace logging = boost::log;
  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;
  
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
      Eigen::Vector3d pos1, pos2, val;
      double volWater=-100, distCurr=-100, distCrit=-100;
    
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
            //std::cerr<<"pos1 " << pos1<<std::endl;
          } else if (i==_cfg->cPos2()) {
            linestream >> pos2[0];
            linestream >> pos2[1];
            linestream >> pos2[2];
            i+=2;
            //std::cerr<<"pos2 " << pos2<<std::endl;
          } else if (i==_cfg->cForc()) {
            linestream >> val[0];
            linestream >> val[1];
            linestream >> val[2];
            i+=2;
            //std::cerr<<"val " << val<<std::endl;
          } else if (i==_cfg->cVolWater()) {
            linestream >> volWater;
            //std::cerr<<"VolWater " << volWater<<std::endl;
          } else if (i==_cfg->cDistCurr()) {
            linestream >> distCurr;
            //std::cerr<<"DistCurr " << distCurr<<std::endl;
          } else if (i==_cfg->cDistCrit()) {
            linestream >> distCrit;
            //std::cerr<<"DistCrit " << distCrit<<std::endl;
          } else {
            linestream >> valD;
          }
        };
        
        if (
             (_particleAll[forceRowNumbTMP]->particleReal(pid1) and 
             _particleAll[forceRowNumbTMP]->particleReal(pid2)) 
             and
             (
               ((_cfg->tF()>=0) and 
                (
                 _particleAll[forceRowNumbTMP]->getP(pid1)->type()==_cfg->tF() or
                 _particleAll[forceRowNumbTMP]->getP(pid2)->type()==_cfg->tF()
                )
               )
               or 
               (_cfg->tF()<0)
             )
           ) {
          std::shared_ptr<force> tmpForce ( new force (_particleAll[forceRowNumbTMP]->getP(pid1), _particleAll[forceRowNumbTMP]->getP(pid2), forceRowNumbTMP, pos1, pos2, val, volWater, distCurr, distCrit));
          _forceRow[forceRowNumbTMP]->addF(tmpForce);
          loadSnap->addForce(tmpForce);
        }
  
      } else if (curLine == _cfg->fAt()) {
        linestream >> valInt;
         BOOST_LOG_SEV(lg, info)<<"Expected forces "<<valInt;
      } else if (curLine == _cfg->nFSt()) {
        linestream >> valInt;
        unsigned int valIntUn = valInt;
        if (valIntUn != loadSnap->timeStep()) {
          BOOST_LOG_SEV(lg, fatal) << "Timestep of force and particle files is not the same!\n"; 
          exit (EXIT_FAILURE);
        }
      }
      curLine++;
    };
    BOOST_LOG_SEV(lg, info)<<_forceRow[forceRowNumbTMP]->elementsNum()<<" forces added.";
    _forceNum+=_forceRow[forceRowNumbTMP]->elementsNum();
};
