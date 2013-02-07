#include "rheometer.h"

rheometer::rheometer(boost::shared_ptr<configopt> cfg, string particlesFileName, string forcesFileName) {
  _cfg = cfg;
  _particlesFileName = particlesFileName;
  _forcesFileName = forcesFileName;
  _particleNum = -1;
  _forceNum = -1;
  
  loadParticles();
  loadForces();
  
  //Create bands
  boost::shared_ptr <bandRow> bandRowTMP (new bandRow(_cfg, _particleAll,  _forceRow));
  boost::shared_ptr <bandRow> _bandRow = bandRowTMP;
  
  
  boost::shared_ptr <exportclass> exp (new exportclass(_cfg, _particleAll));
  exp->exportVTK();
};

void rheometer::loadParticles() {
  std::ifstream _file;
  _file.open(_particlesFileName.c_str());
  
  std::string   line;
  int curLine = 1;
  unsigned long long maxId = 0;
  std::vector <boost::shared_ptr<particle> > tmpPartVector;

  while(std::getline(_file, line)) {
    std::stringstream linestream(line);
    std::string data;
    
    int valInt;
    double valD;
    double pR;
    int pT;
    unsigned long long pId;
    Eigen::Vector3d pC, pV, pO;
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
        } else {
          linestream >> valD;
        }
      }
      maxId = max(pId, maxId);
      boost::shared_ptr<particle> tmpParticle ( new particle (pId, pT, pR, pC,pV, pO));
      tmpPartVector.push_back(tmpParticle);

    } else if (curLine == _cfg->nAt()) {
      linestream >> valInt;
      _particleNum = valInt;
      std::cerr<<"Expected number of particles "<<_particleNum<<std::endl;
    }
    curLine++;
  };
  boost::shared_ptr<particleRow> particleTMP ( new particleRow(maxId+1));
  _particleAll = particleTMP;
  
  for(std::vector<boost::shared_ptr<particle> >::iterator it = tmpPartVector.begin(); it != tmpPartVector.end(); ++it) {
    _particleAll->addP(*it);
  }
  std::cerr<<_particleAll->elementsNum()<<" particles added"<<std::endl;
};

void rheometer::loadForces() {
  std::ifstream _file;
  _file.open(_forcesFileName.c_str());
  
  std::string   line;
  int curLine = 1;
  
  //std::vector <boost::shared_ptr<particle> > tmpPartVector;
  
  int valInt;
  double valD;
  
  boost::shared_ptr<forceRow> forceTMP ( new forceRow());
  _forceRow = forceTMP;
  
  
  
  while(std::getline(_file, line)) {
    std::stringstream linestream(line);
    std::string data;
    
    unsigned long long pid1, pid2;
    Eigen::Vector3d pos1, pos2, val;
  
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
      
      
      boost::shared_ptr<force> tmpForce ( new force (pid1, pid2, pos1, pos2, val));
      _forceRow->addF(tmpForce);

    } else if (curLine == _cfg->fAt()) {
      linestream >> valInt;
      _forceNum = valInt;
      std::cerr<<"Expected number of forces "<<_forceNum<<std::endl;
    }
    curLine++;
  };
  
  std::cerr<<_forceRow->elementsNum()<<" forces added"<<std::endl;
  
};
