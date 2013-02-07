#ifndef RHEOMETERCLASS
#define RHEOMETERCLASS

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "config.h"
#include "particle.h"
#include "band.h"
#include "force.h"
#include "export.h"

using namespace std;
class rheometer {
  private:
    boost::shared_ptr<configopt> _cfg;
    string _particlesFileName;
    string _forcesFileName;
    
    unsigned long long  _particleNum;
    unsigned long long _forceNum;
    boost::shared_ptr <particleRow> _particleAll;
    boost::shared_ptr <forceRow> _forceRow;
    boost::shared_ptr <bandRow> _bandRow;
  public:
    rheometer(boost::shared_ptr<configopt>, string, string);
    void loadParticles();
    void loadForces();
};

#endif
