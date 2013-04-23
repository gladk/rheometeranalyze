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
#include "math_custom.h"
#include <iostream>

particle::particle(unsigned long long id, int type, unsigned int fileid, double rad, double mass, double dens, Eigen::Vector3d c, Eigen::Vector3d v, Eigen::Vector3d o) {
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
  _stressTensor = _stressTensor.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
};

particle::particle() {
  _id = -1;
  _type = -1;
  _rad = -1;
  _d = -1;
  _m = -1;
  _fileId = -1;
  _c = Eigen::Vector3d::Zero();
  _v = Eigen::Vector3d::Zero();
  _o = Eigen::Vector3d::Zero();
  _vZylindrical = Eigen::Vector3d::Zero();
  _disable = false;
  _axisMatrix = _axisMatrix.Zero();
  _stressTensor = _stressTensor.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
};

void particle::addStress(Eigen::Matrix3d addStressTensor) {
  _stressTensor += addStressTensor;
};

Eigen::Matrix3d particle::stressTensor() {
  return _stressTensor;
};

Eigen::Matrix3d particle::stressTensorAVG() {
  if (_contactParticles.size()>0) {
    return _stressTensor/this->vol();
  } else {
    return _stressTensor.Zero();
  }
};

double particle::realAngular() {
  if (not(_calculateVel)) { 
    std::cerr << "Zylindrical velocities are not set yet! Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  };
  return _vZylindrical(2)/_posZyl(0);
};

double particle::stressPress() {
  Eigen::Matrix3d stressTMP = this->kinEnergie() + this->stressTensor();
  return stressTMP.trace()/3.0;
};

double particle::stressTau() {
  Eigen::Matrix3d stressTMP = this->stressTensor();
  return sqrt(stressTMP(1)*stressTMP(1) + stressTMP(2)*stressTMP(2) + stressTMP(5)*stressTMP(5));;
};

Eigen::Matrix3d particle::kinEnergie() {
  if (not(_calculateVel)) { 
    std::cerr << "Zylindrical velocities are not set yet! Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  };
  return _m*_vZylindrical*_vZylindrical.transpose();
};

Eigen::Vector3d particle::vZyl() {
  if (not(_calculateVel)) { 
    std::cerr << "Zylindrical velocities are not set yet! Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  };
  return _vZylindrical;
};

unsigned int particle::contacts() {
  return _contactParticles.size();
};

void particle::setPosZyl(Eigen::Vector3d zyl, Eigen::Quaternion<double> rotateCCz) {
  _posZyl = zyl;
  
  double const& rho = zyl(0);
  double const& z = zyl(1);
  double const& phi = zyl(2);
  
  Eigen::Vector3d dr = Eigen::Vector3d(cos(phi), sin(phi), 0.0); dr = rotateCCz*dr;
  Eigen::Vector3d df = Eigen::Vector3d(-sin(phi), cos(phi), 0.0); df = rotateCCz*df;
  Eigen::Vector3d dz = Eigen::Vector3d(0.0, 0.0, z); dz = rotateCCz*dz;
  
  dr.normalize(); dz.normalize(); df.normalize();
  _axisMatrix << dr, dz, df;
  _axisMatrix.transposeInPlace();
  
  Eigen::Matrix3d valTempMatrix; valTempMatrix << _v, _v, _v;
  valTempMatrix.transposeInPlace();
  valTempMatrix = _axisMatrix.cwiseProduct(valTempMatrix);
  
  _vZylindrical = Eigen::Vector3d(valTempMatrix.row(0).sum(),
                                  valTempMatrix.row(1).sum(),
                                  valTempMatrix.row(2).sum());
  _calculateVel = true;
};

void particle::addParticleContact(std::shared_ptr<particle> addParticle) {
  if (_id==addParticle->id()) {
    std::cerr << "Cannot add the contact force between the same particles "<<_id<< " and "<< addParticle->id() << " has been added 2 times! Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  }
  for (unsigned int i = 0; i < _contactParticles.size(); i++) {
    if (_contactParticles[i]->id()==addParticle->id()) {
      std::cerr << "The Force between particles "<<_contactParticles[i]->id()<< " and "<< addParticle->id() << " has been added 2 times! Exiting..."<<std::endl; 
      exit (EXIT_FAILURE);
    }
  }
  _contactParticles.push_back(addParticle);
};
