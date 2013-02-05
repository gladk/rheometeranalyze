#include "rheometer.h"

rheometer::rheometer(boost::shared_ptr<configopt> cfg, string particlesFileName) {
  _cfg = cfg;
  _particlesFileName = particlesFileName;
  loadParticles();
  boost::shared_ptr<particleRow> part_s (new particleRow(5));

};

void rheometer::loadParticles() {
  std::ifstream _file;
  _file.open(_particlesFileName.c_str());
  
  std::string   line;
  while(std::getline(_file, line)) {
    std::stringstream linestream(line);
    std::string data;
    double valX;
    double valY;
    
    for (int i=0; i++; i<_cfg->maxC()) {
      linestream >> valX;
    }
  }

};
