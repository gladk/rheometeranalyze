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

#include "rheometer.h"
#include "snapshot.h"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/algorithm/string.hpp>

rheometer::rheometer(std::shared_ptr<configopt> cfg) {
  _cfg = cfg;
  _particleNum = 0;
  _forceNum = 0;
  _snapshots = _cfg->snapshot();
  std::shared_ptr <exportclass> exp;
  
  //Create bands
  if (_cfg->increment()) {
    std::vector<std::shared_ptr<bandRow>> allBandRows;
    for(unsigned long long i=0; i<_snapshots->size(); i++) {
      loadParticles(i);
      _bandRow = std::make_shared<bandRow>(_cfg, _particleAll,  _forceRow);
      calculateLocalDeformations();
      
      _bandRow->clear();
      _particleAll.clear(); // _particleAll.shrink_to_fit();
      _forceRow.clear(); //_forceRow.shrink_to_fit();
      _snapshots->clear();
      allBandRows.push_back(_bandRow);
    }
    _bandRow = std::make_shared<bandRow>(_cfg, allBandRows);
    exp = std::make_shared<exportclass>(_cfg, _bandRow);
  } else {
    loadParticles();
    _bandRow = std::make_shared<bandRow>(_cfg, _particleAll,  _forceRow);
     calculateLocalDeformations();
    
    exp = std::make_shared<exportclass>(_cfg, _bandRow, _forceRow);
  }
  
  // ==========================ForceChain
  
  exp->forceChain();
  
  // ==========================ForceChain
  
  if (_cfg->Vtk()) exp->VTK();
  if (_cfg->Utwente()) exp->Utwente();
  
  exp->gnuplotSchearRate();
  
  exp->torque();
  
  if (_cfg->contact()) {
    exp->gnuplotContactAnalyze(100);
    exp->gnuplotContactNumberAnalyze();
    exp->gnuplotContactWet();
  }
  
  if (_cfg->followContact()) {
    exp->gnuplotContactFollow();
  }
  
  if (_cfg->intOri()>0) {
    exp->intOri();
  }
  
  if (_cfg->wetParticle()>0) {
    exp->gnuplotWetParticles(_cfg->wetParticle());
  }
  
};

void rheometer::loadParticles(unsigned long long numLoad) {
  namespace src = boost::log::sources;
  namespace logging = boost::log;
  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;
    
  unsigned int partNumbCounter  = 1;
  
  unsigned long long loadSnapshotsNum = _snapshots->size();
  
  if (_cfg->increment()) {
    std::cerr<<"numLoad: "<<numLoad<<std::endl;
    loadSnapshotsNum=numLoad+1;
  }
  
  for(unsigned int i=numLoad; i<loadSnapshotsNum; i++) {
    
    std::shared_ptr<snapshot> snapshotCur = _snapshots->getSnapshot(i);
    
    std::ifstream _file;
    boost::iostreams::filtering_istream in;
    _file.open(snapshotCur->getParticleFile().string(), std::ios_base::in | std::ios_base::binary);
    
    if ((snapshotCur->getParticleFile()).extension().string()==".bz2") {
      in.push(boost::iostreams::bzip2_decompressor());
    } else if ((snapshotCur->getParticleFile()).extension().string()==".gz") {
      in.push(boost::iostreams::gzip_decompressor());
    }
    
    in.push(_file);
    std::string  line;
    int curLine = 1;
    unsigned long long maxId = 0;
    unsigned long long expectedParticles = 0;
    std::vector <std::shared_ptr<particle> > tmpPartVector;
    
    //======================================================================
    //======================================================================
    while(std::getline(in, line)) {
      boost::algorithm::trim(line);
      if (not(line.empty())) {
        std::stringstream linestream(line);
        std::string data;
        
        int valInt=0;
        double valD=0;
        double pR, pM, pD, volWater=0.0;
        int pT=0;
        unsigned long long pId=0;
        Eigen::Vector3d pC=Eigen::Vector3d::Zero();
        Eigen::Vector3d pV=Eigen::Vector3d::Zero();
        Eigen::Vector3d pO=Eigen::Vector3d::Zero();
        
        bool addId = false;   // flag shows, that Id has been added
        bool addC  = false;   // flag shows, that center coordinates have been added
        
        if (curLine>=_cfg->nDat()) {
          for (int i=1; i<=_cfg->maxC(); i++) {
            if (i==_cfg->cId()) {
              linestream >> pId;
              addId = true;
            } else if (i==_cfg->cT()) {
              linestream >> pT;
              //std::cerr<<pT;
            } else if (i==_cfg->cC()) {
              linestream >> pC[0];
              linestream >> pC[1];
              linestream >> pC[2];
              i+=2;
              addC = true;
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
          
          if ((((_cfg->tC()>=0) and (pT == _cfg->tC())) or (_cfg->tC()<0)) and 
             ((expectedParticles>0) and (tmpPartVector.size()<expectedParticles)) and
             (addId and addC)) {
            maxId = max(pId, maxId);
            std::shared_ptr<particle> tmpParticle = std::make_shared<particle>(pId, pT, partNumbCounter-1, pR, pM, pD, pC,pV, pO, volWater);
            if (_cfg->intOri() > 0) tmpParticle->createIntOri(_cfg->intOri());
            tmpPartVector.push_back(tmpParticle);
            snapshotCur->addParticle(tmpParticle);
          }
    
        } else if (curLine == _cfg->nAt()) {
          linestream >> expectedParticles;
          unsigned int FileNumber = partNumbCounter;
          if (_cfg->increment()) {FileNumber = numLoad+1;}
          BOOST_LOG_SEV(lg, info)<<"File "<<FileNumber<<"/"<<_snapshots->size()<< " (" << snapshotCur->getParticleFile() << "); " <<"Expected particles "<<expectedParticles;
        } else if (curLine == _cfg->nPSt()) {
          linestream >> valInt;
          snapshotCur->setTimeStep(valInt, _cfg->dT());
        }
      }
      curLine++;
    };
    //======================================================================
    //======================================================================
    
    std::shared_ptr<particleRow> particleTMP = std::make_shared<particleRow>(maxId+1);
    _particleAll.push_back(particleTMP);
    unsigned int partNumbTMP  = _particleAll.size()-1;
    
    for(auto p : tmpPartVector) {
       _particleAll[partNumbTMP]->addP(p);
    }
    BOOST_LOG_SEV(lg, info)<<_particleAll[partNumbTMP]->elementsNum()<<" particles added.";
    _particleNum+=_particleAll[partNumbTMP]->elementsNum();
    partNumbCounter++;
    
    this->loadForces(snapshotCur);
  }
  
  if (!_cfg->increment()) {
    _snapshots->sortRow();
  }
  BOOST_LOG_SEV(lg, info)<<"The total number of added particles is "<<_particleNum;
  BOOST_LOG_SEV(lg, info)<<"The total number of added forces is "<<_forceNum;
};


void rheometer::loadForces(std::shared_ptr<snapshot> loadSnap) {
  namespace src = boost::log::sources;
  namespace logging = boost::log;
  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;
  
    std::ifstream _file;
    boost::iostreams::filtering_istream in;
    _file.open(loadSnap->getForceFile().string(), std::ios_base::in | std::ios_base::binary);
    
    if ((loadSnap->getForceFile()).extension().string()==".bz2") {
      in.push(boost::iostreams::bzip2_decompressor());
    } else if ((loadSnap->getForceFile()).extension().string()==".gz") {
      in.push(boost::iostreams::gzip_decompressor());
    }
    
    in.push(_file);
    std::string   line;
    int curLine = 1;
    
    int valInt;
    double valD;
    
    std::shared_ptr<forceRow> forceTMP = std::make_shared<forceRow>();
    _forceRow.push_back(forceTMP);
    unsigned int forceRowNumbTMP  = _forceRow.size()-1;
    
    
    while(std::getline(in, line)) {
      boost::algorithm::trim(line);
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
          std::shared_ptr<force> tmpForce = std::make_shared<force>(_particleAll[forceRowNumbTMP]->getP(pid1), _particleAll[forceRowNumbTMP]->getP(pid2), forceRowNumbTMP, pos1, pos2, val, volWater, distCurr, distCrit);
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

void rheometer::calculateLocalDeformations() {
  double dT;
  if (_cfg->_timeCur < 0) {
    dT = (_snapshots->timeMax() - _snapshots->timeMin())/2.0;
    _cfg->_timeCur = dT;
  } else {
    dT = (_snapshots->timeMax() - _snapshots->timeMin());
    _cfg->_timeCur += dT;
  }
  _cfg->_timeStepCur = _snapshots->timeStepAvg();
  
  for(unsigned int b=0; b<_bandRow->size(); b++) {
    std::shared_ptr<band> bandTMP = _bandRow->getBand(b);
    
    /* Calculate local strain, deformations according to Sakaie(5)
     * const double r = _cfg->Dout()/2.0;
     * const double N = _bandRow->omega0AVG()*_snapshots->timeAvg()/(2*M_PI);
     * const double gamma =  TwoPiNr*bandTMP->dOmegadR();
    */
    
    const double gamma =  _cfg->_gamma[b] + bandTMP->scherRate()*dT;
    bandTMP->gamma(gamma);
    _cfg->_gamma[b] = gamma;
  }
}
