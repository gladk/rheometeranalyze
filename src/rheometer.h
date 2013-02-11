#ifndef RHEOMETERCLASS
#define RHEOMETERCLASS

#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include <Eigen/Dense>
#include "config.h"
#include "particle.h"
#include "band.h"
#include "force.h"
#include "export.h"

using namespace std;
class rheometer {
  private:
    std::shared_ptr<configopt> _cfg;
    string _particlesFileName;
    string _forcesFileName;
    
    unsigned long long  _particleNum;
    unsigned long long _forceNum;
    std::shared_ptr <particleRow> _particleAll;
    std::shared_ptr <forceRow> _forceRow;
    std::shared_ptr <bandRow> _bandRow;
  public:
    rheometer(std::shared_ptr<configopt>, string, string);
    void loadParticles();
    void loadForces();
};

#endif
