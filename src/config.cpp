/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014 Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RheometerAnalyze is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RheometerAnalyze.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "config.h"

namespace fs = boost::filesystem;

configopt::configopt(const std::string &filename) {
  boost::property_tree::ptree pt;
  boost::property_tree::ini_parser::read_ini(filename, pt);
  _maxC = -1;
  _maxCF = -1;
  _wetParticle = -1;
  
  _c = Eigen::Vector3d(pt.get<double>("Rheometer.cX"), pt.get<double>("Rheometer.cY"), pt.get<double>("Rheometer.cZ"));
  _o = Eigen::Vector3d(pt.get<double>("Rheometer.oX"), pt.get<double>("Rheometer.oY"), pt.get<double>("Rheometer.oZ"));
  _o.normalize();

  _Din = pt.get<double>("Rheometer.Din");
  _Dout = pt.get<double>("Rheometer.Dout");
  _H = pt.get<double>("Rheometer.H");

  _SecRadial = pt.get<int>("Rheometer.SecRadial");
  _SecZ = pt.get<int>("Rheometer.SecZ");
  _SecFi = pt.get<int>("Rheometer.SecFi");
  
  _dT = pt.get<double>("Rheometer.dT");
  
  _g = Eigen::Vector3d(pt.get<double>("Rheometer.gX"), pt.get<double>("Rheometer.gY"), pt.get<double>("Rheometer.gZ"));
  _g.normalize();
  
  _aS = pt.get<double>("Rheometer.aS"); _aS *=M_PI/180.0;
  _aE = pt.get<double>("Rheometer.aE"); _aE *=M_PI/180.0;
  if (_aS>_aE) {
    std::cerr<<"aS cannot be more, than aE. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }

  _nAt = pt.get<int>("Particle.nAt");
  _nPSt = pt.get<int>("Particle.nPSt");
  _nDat = pt.get<int>("Particle.nDat");
  _tC = pt.get<double>("Particle.tC");
  _tR = pt.get<double>("Particle.tR");

  _cId = pt.get<int>("Particle.cId"); _maxColumnCheck(_cId, 0);
  _cT = pt.get<int>("Particle.cT"); _maxColumnCheck(_cT, 0);
  _cC = pt.get<int>("Particle.cC"); _maxColumnCheck(_cC, 2);
  _cV = pt.get<int>("Particle.cV"); _maxColumnCheck(_cV, 2);
  _cO = pt.get<int>("Particle.cO"); _maxColumnCheck(_cO, 2);
  _cR = pt.get<double>("Particle.cR"); _maxColumnCheck(_cR, 0);
  _cM = pt.get<double>("Particle.cM"); _maxColumnCheck(_cM, 0);
  _cD = pt.get<double>("Particle.cD"); _maxColumnCheck(_cD, 0);
  _cVolWaterP = pt.get<double>("Particle.cVolWaterP"); _maxColumnCheck(_cVolWaterP, 0);
  
  
  //Force
  _fAt = pt.get<int>("Force.nFt");
  _nFSt = pt.get<int>("Force.nFSt");
  _fDat = pt.get<int>("Force.nDat");
  _cPos1 = pt.get<int>("Force.cPos1"); _maxColumnCheckForce(_cPos1, 2);
  _cPos2 = pt.get<int>("Force.cPos2"); _maxColumnCheckForce(_cPos2, 2);
  _cPos1ID = pt.get<int>("Force.cPos1ID"); _maxColumnCheckForce(_cPos1ID, 0);
  _cPos2ID = pt.get<int>("Force.cPos2ID"); _maxColumnCheckForce(_cPos2ID, 0);
  _cForc = pt.get<int>("Force.cForc"); _maxColumnCheckForce(_cForc, 2);
  _cVolWater = pt.get<int>("Force.cVolWater"); _maxColumnCheckForce(_cVolWater, 0);
  _cDistCurr = pt.get<int>("Force.cDistCurr"); _maxColumnCheckForce(_cDistCurr, 0);
  _cDistCrit = pt.get<int>("Force.cDistCrit"); _maxColumnCheckForce(_cDistCrit, 0);
  _tF = pt.get<int>("Force.tF");
  
  _utwente = false;
  _vtk = 0;
  _contact = false;
  _followContact = false;
  _intOri = -1;
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
