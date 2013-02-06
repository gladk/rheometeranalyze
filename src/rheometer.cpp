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
  while(std::getline(_file, line)) {
    std::stringstream linestream(line);
    std::string data;
    
    int valInt;
    double valD;

    double pR;
    int pId, pT;
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
    } else if (curLine == _cfg->nAt()) {
      linestream >> valInt;
      _particleNum = valInt;
    }
    curLine++;
  }

};
