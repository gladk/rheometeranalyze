/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "particle.h"
#include <iostream>

particle::particle(unsigned long long id, int type, unsigned int fileid, double rad, double mass, double dens, Eigen::Vector3f c, Eigen::Vector3f v, Eigen::Vector3f o) {
  _id = id;
  _type = type;
  _rad = rad;
  _fileId = fileid;
  _c = c;
  _v = v;
  _o = o;
  _m = mass;
  _d = dens;
  _disable = false;
  _axisMatrix = _axisMatrix.Zero();
  _velMatrix = _velMatrix.Zero();
  _stressTensor = _stressTensor.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
  _contacts = 0;
};

particle::particle() {
  _id = -1;
  _type = -1;
  _rad = -1;
  _d = -1;
  _m = -1;
  _fileId = -1;
  _c = Eigen::Vector3f::Zero();
  _v = Eigen::Vector3f::Zero();
  _o = Eigen::Vector3f::Zero();
  _vZylindrical = Eigen::Vector3f::Zero();
  _disable = false;
  _axisMatrix = _axisMatrix.Zero();
  _velMatrix = _velMatrix.Zero();
  _stressTensor = _stressTensor.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
  _contacts = 0;
};

void particle::set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df) {
  _axisMatrix = _axisMatrix.Zero();
  dr.normalize(); dz.normalize(); df.normalize();
  _axisMatrix << dr, dz, df;
  _axisMatrix.transposeInPlace();
  if (not(_calculateVel)) { calculateVel();}
};

void particle::calculateVel() {
  Eigen::Matrix3f velTempMatrix; velTempMatrix << _v, _v, _v;
  _vZylindrical = Eigen::Vector3f(_v.dot(this->dr()), _v.dot(this->dz()), _v.dot(this->df()));
  velTempMatrix.transposeInPlace();
  _velMatrix = _axisMatrix.cwiseProduct(velTempMatrix);
  _calculateVel = true;
};

void particle::addStress(Eigen::Matrix3f addStressTensor) {
  _contacts +=1;
  _stressTensor += addStressTensor;
};

Eigen::Matrix3f particle::stressTensor() {
  return _stressTensor;
};

Eigen::Matrix3f particle::stressTensorAVG() {
  if (_contacts>0) {
    return _stressTensor/this->vol();
  } else {
    return _stressTensor.Zero();
  }
};

double particle::realAngular() {
  if (not(_calculateVel)) { calculateVel();};
  return _vZylindrical(2)/_posZyl(0);
};

double particle::stressPress() {
  Eigen::Matrix3f stressTMP = this->kinEnergie() + this->stressTensor();
  return stressTMP.trace()/3.0;
};

double particle::stressTau() {
  Eigen::Matrix3f stressTMP = this->stressTensor();
  return sqrt(stressTMP(1)*stressTMP(1) + stressTMP(2)*stressTMP(2) + stressTMP(5)*stressTMP(5));;
};


Eigen::Matrix3f particle::kinEnergie() {
  if (not(_calculateVel)) { calculateVel();};
  return _m*_vZylindrical*_vZylindrical.transpose();
};

Eigen::Vector3f particle::vZyl() {
  if (not(_calculateVel)) { calculateVel();};
  return _vZylindrical;
};
