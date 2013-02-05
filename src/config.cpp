#include "config.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include <iostream>

configopt::configopt(const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  std::cout << pt.get<double>("Rheometer.cX") << std::endl;
  std::cout << pt.get<double>("Rheometer.cY") << std::endl;
  std::cout << pt.get<double>("Rheometer.cZ") << std::endl;
  std::cout << pt.get<double>("Rheometer.Din") << std::endl;
  std::cout << pt.get<double>("Rheometer.Dout") << std::endl;
  std::cout << pt.get<int>("Rheometer.SecRadial") << std::endl;
  std::cout << pt.get<int>("Rheometer.SecZ") << std::endl;
}
