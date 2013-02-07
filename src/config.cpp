#include "config.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

configopt::configopt(const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  _maxC = -1;
  _maxCF = -1;
  
  _c = Eigen::Vector3d(pt.get<double>("Rheometer.cX"), pt.get<double>("Rheometer.cY"), pt.get<double>("Rheometer.cZ"));
  _o = Eigen::Vector3d(pt.get<double>("Rheometer.oX"), pt.get<double>("Rheometer.oY"), pt.get<double>("Rheometer.oZ"));
  _o.normalize();

  _Din = pt.get<double>("Rheometer.Din");
  _Dout = pt.get<double>("Rheometer.Dout");
  _H = pt.get<double>("Rheometer.H");

  _SecRadial = pt.get<int>("Rheometer.SecRadial");
  _SecZ = pt.get<int>("Rheometer.SecZ");
  
  _g = Eigen::Vector3d(pt.get<double>("Rheometer.gX"), pt.get<double>("Rheometer.gY"), pt.get<double>("Rheometer.gZ"));
  _g.normalize();

  _nAt = pt.get<int>("Particle.nAt");
  _nDat = pt.get<int>("Particle.nDat");

  _cId = pt.get<int>("Particle.cId"); _maxColumnCheck(_cId, 0);
  _cT = pt.get<int>("Particle.cT"); _maxColumnCheck(_cT, 0);
  _cC = pt.get<int>("Particle.cC"); _maxColumnCheck(_cC, 2);
  _cV = pt.get<int>("Particle.cV"); _maxColumnCheck(_cV, 2);
  _cO = pt.get<int>("Particle.cO"); _maxColumnCheck(_cO, 2);
  _cR = pt.get<double>("Particle.cR"); _maxColumnCheck(_cR, 0);
  
  
  //Force
  _fAt = pt.get<int>("Force.nFt");
  _fDat = pt.get<int>("Force.nDat");
  _cPos1 = pt.get<int>("Force.cPos1"); _maxColumnCheckForce(_cPos1, 2);
  _cPos2 = pt.get<int>("Force.cPos2"); _maxColumnCheckForce(_cPos2, 2);
  _cPos1ID = pt.get<int>("Force.cPos1ID"); _maxColumnCheckForce(_cPos1ID, 0);
  _cPos2ID = pt.get<int>("Force.cPos2ID"); _maxColumnCheckForce(_cPos2ID, 0);
  _cForc = pt.get<int>("Force.cForc"); _maxColumnCheckForce(_cForc, 2);
  
};

void configopt::_maxColumnCheck(int col, int addN) {
  if (col+addN > _maxC) {
    _maxC = col+addN;
  }
};

void configopt::_maxColumnCheckForce(int col, int addN) {
  if (col+addN > _maxCF) {
    _maxCF = col+addN;
  }
};
