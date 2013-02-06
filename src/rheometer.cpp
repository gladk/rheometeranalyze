#include "rheometer.h"

rheometer::rheometer(boost::shared_ptr<configopt> cfg, string particlesFileName) {
  _cfg = cfg;
  _particlesFileName = particlesFileName;
  _particleNum = -1;
  loadParticles();
  boost::shared_ptr<particleRow> part_s (new particleRow(5));

};

void rheometer::loadParticles() {
  std::ifstream _file;
  _file.open(_particlesFileName.c_str());
  
  std::string   line;
  int curLine = 1;
  long long maxId = -1;
  std::vector <boost::shared_ptr<particle> > tmpPartVector;

  while(std::getline(_file, line)) {
    std::stringstream linestream(line);
    std::string data;
    
    int valInt;
    double valD;
    double pR;
    int pT;
    long long pId;
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
    }
    curLine++;
  };
  boost::shared_ptr<particleRow> particleTMP ( new particleRow(maxId+1));
  _particleAll = particleTMP;
  
  for(std::vector<boost::shared_ptr<particle> >::iterator it = tmpPartVector.begin(); it != tmpPartVector.end(); ++it) {
    _particleAll->addP(*it);
  }
  std::cerr<<_particleAll->elementsNum()<<" particles added"<<std::endl;
  
  //Create bands
  boost::shared_ptr <bandRow> bandRowTMP (new bandRow(_cfg, _particleAll));
  boost::shared_ptr <bandRow> _bandRow = bandRowTMP;
  

};
