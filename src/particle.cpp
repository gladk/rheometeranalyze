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

#include "particle.h"
#include "math_custom.h"
#include <iostream>
#include <boost/foreach.hpp>

particle::particle(unsigned long long id, int type, unsigned int fileid, 
                   double rad, double mass, double dens, Eigen::Vector3d c, 
                   Eigen::Vector3d v, Eigen::Vector3d o, double volwater) {
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
  _stressTensorCap = _stressTensorCap.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
  _sizeIntOri = -1;
  _shearBand = false;
  _highStressed = -1;
  _snapshot = -1;
  _volwater = volwater;
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
  _stressTensorCap = _stressTensorCap.Zero();
  _posZyl = _posZyl.Zero();
  _calculateVel = false;
  _press = 0.0;
  _tau = 0.0;
  _sizeIntOri = -1;
  _shearBand = false;
  _highStressed = -1;
  _snapshot = -1;
};

void particle::addStress(Eigen::Matrix3d addStressTensor) {
  _stressTensor += addStressTensor;
};

void particle::addStressCap(Eigen::Matrix3d addStressTensor) {
  _stressTensorCap += addStressTensor;
};

Eigen::Matrix3d particle::stressTensor() {
  return _stressTensor;
};

Eigen::Matrix3d particle::stressTensorCap() {
  return _stressTensorCap;
};

Eigen::Matrix3d particle::stressTensorAVG() {
  if (_contactParticles.size()>0) {
    return _stressTensor/this->vol();
  } else {
    return _stressTensor.Zero();
  }
};

Eigen::Matrix3d particle::stressTensorCapAVG() {
  if (_contactParticlesWet.size()>0) {
    return _stressTensorCap/this->vol();
  } else {
    return _stressTensorCap.Zero();
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
  /*
   * 
   * The formula (15, part1) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
   * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
   * 
   */
  return _m*_vZylindrical*_vZylindrical.transpose();
};

double particle::kinEnergieDouble() {
  return _m*_v.norm()*_v.norm()*0.5;
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

unsigned int particle::wetContacts() {
  return _contactParticlesWet.size();
};

double particle::wetContactsAverageDistance() {
  double wetContactsAverageDistance = 0.0;
  for (unsigned int i = 0; i < _contactParticlesWet.size(); i++) {
    std::shared_ptr<particle> tmpParticle = _contactParticles[i];
    double tmpDist = (_c - tmpParticle->c()).norm() - (_rad + tmpParticle->rad());
    wetContactsAverageDistance += tmpDist;
  }
  if (_contactParticlesWet.size()>0) {
    return wetContactsAverageDistance/_contactParticlesWet.size();
  } else {
    return 0.0;
  }
};

void particle::setLocalCoord(Eigen::Vector3d loc, Eigen::Quaternion<double> rotateCC)  {
  _posZyl = cart_to_cyl(rotateCC*loc);
  _axisMatrix = get_axes_coord(_posZyl, rotateCC);
  
  _vZylindrical = get_cyl_rotated_vector(_v, _axisMatrix);
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
  if (_sizeIntOri > 0) {
    Eigen::Vector3d tmpL = addParticle->c()-_c;
    Eigen::Vector3d lZ = get_cyl_rotated_vector(tmpL, _axisMatrix);
    _normContOri->addinteraction(lZ);
  }
};

void particle::addParticleContactWet(std::shared_ptr<particle> addParticle) {
  this->addParticleContact(addParticle);
  _contactParticlesWet.push_back(addParticle);
  if (_sizeIntOri > 0) {
    Eigen::Vector3d tmpL = addParticle->c()-_c;
    Eigen::Vector3d lZ = get_cyl_rotated_vector(tmpL, _axisMatrix);
    _capiContOri->addinteraction(lZ);
  }
};

void particle::createIntOri(unsigned short intNumb) {
  _sizeIntOri = intNumb;
  _normContOri = std::make_shared<interactionori>(intNumb);
  _capiContOri = std::make_shared<interactionori>(intNumb);;
};

double particle::stressSigma1() {
  const double sigma1 = ((_stressTensor(0) + _stressTensor(8))/2.0 + 
    sqrt(pow((_stressTensor(0) - _stressTensor(8))/2.0, 2) + 
    _stressTensor(2)*_stressTensor(2)))/this->vol();
  return sigma1;
};

double particle::stressSigma3() {
  const double sigma3 = ((_stressTensor(0) + _stressTensor(8))/2.0 - 
    sqrt(pow((_stressTensor(0) - _stressTensor(8))/2.0, 2) + 
    _stressTensor(2)*_stressTensor(2)))/this->vol();
  return sigma3;
};

bool particle::wetContact(std::shared_ptr <particle> p) {
  BOOST_FOREACH(std::shared_ptr <particle> j,  _contactParticlesWet) {
    if (p==j) return true;
  }
  return false;
};

int particle::highStressedContacts() {
  int highStressedContactsTMP = 0;
  BOOST_FOREACH(std::shared_ptr <particle> p,  _contactParticles) {
    if (this->highStressedContact(p) and not(p->disabled()) and not(this->disabled())) highStressedContactsTMP++;
  }
  return highStressedContactsTMP;
};

bool particle::highStressedContact(std::shared_ptr <particle> p) {
  BOOST_FOREACH(std::shared_ptr <particle> j,  _contactParticles) {
    if ((p==j) and (p->highStress()>0) and 
        not(this->wetContact(p))) return true;
  }
  return false;
};

double particle::relativeWaterVol() {
  if (this->volwater()>0) {
    return this->volwater()/this->vol();
  } else {
    return 0;
  }
};

InteractionsMatrix particle::normContOri() {
  return _normContOri->interactions();
};

InteractionsMatrix particle::capiContOri() {
  return _capiContOri->interactions();
};
