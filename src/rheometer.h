#ifndef RHEOMETERCLASS
#define RHEOMETERCLASS

#include <string>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "config.h"
#include "particle.h"

using namespace std;
class rheometer {
  private:
    boost::shared_ptr<configopt> _cfg;
    string _particlesFileName;
  public:
    rheometer(boost::shared_ptr<configopt>, string);
    void loadParticles();
};

#endif
