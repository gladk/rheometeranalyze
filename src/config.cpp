#include "config.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

configopt::configopt(const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  
  _c = Eigen::Vector3d(pt.get<double>("Rheometer.cX"), pt.get<double>("Rheometer.cY"), pt.get<double>("Rheometer.cZ"));
  _o = Eigen::Vector3d(pt.get<double>("Rheometer.oX"), pt.get<double>("Rheometer.oY"), pt.get<double>("Rheometer.oZ"));
  _o.normalize();

  _Din = pt.get<double>("Rheometer.Din");
  _Dout = pt.get<double>("Rheometer.Dout");
  _SecRadial = pt.get<int>("Rheometer.SecRadial");
  _SecZ = pt.get<int>("Rheometer.SecZ");

};

