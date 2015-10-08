/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014, 2015 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014, 2015 Anton Gladky <gladky.anton@gmail.com>

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

#include <boost/foreach.hpp>
#include "snapshot.h"

snapshot::snapshot(fs::path particlesFileName, fs::path forcesFileName, 
                   unsigned long long timeStep) {
  _particlesFileName = particlesFileName;
  _forcesFileName = forcesFileName;
  _timeStep = timeStep;
  _id=-1;
};

void snapshot::setParticlesFileName(fs::path particlesFileName) {
  _particlesFileName = particlesFileName;
};

void snapshot::setForcesFileName(fs::path forcesFileName) {
  _forcesFileName = forcesFileName;
};

void snapshot::setTimeStep(const unsigned long long timeStep, const double dT) {
  _timeStep = timeStep;
  _time = _timeStep*dT;
};

unsigned long long snapshot::timeStep() const {
  return _timeStep;
};

unsigned long long snapshot::size() const {
  return _particles.size();
};

double snapshot::time() const {
  return _time;
};

fs::path snapshot::getParticleFile() {
  return _particlesFileName;
};

fs::path snapshot::getForceFile() {
  return _forcesFileName;
};

void snapshot::addParticle(std::shared_ptr<particle> particleTmp) {
  _particles.push_back(particleTmp);
};

void snapshot::addForce(std::shared_ptr<force> forceTmp) {
  _forces.push_back(forceTmp);
};

void snapshot::id(unsigned short id) {
  _id = id;
  BOOST_FOREACH(std::shared_ptr<particle> p, _particles) {
    p->snapshot(id);
  }
};

unsigned short snapshot::id() const {
  return _id;
};

std::vector <std::shared_ptr<particle> >  snapshot::particles() {
  return _particles;
};

std::vector <std::shared_ptr<force> >  snapshot::forces() {
  return _forces;
};

double snapshot::torque(Eigen::Vector3d rotationAxis, Eigen::Vector3d zeroPoint, int typeAnalyze) {
  double torqueRet = 0.0;
  rotationAxis.normalize();
  
  BOOST_FOREACH(std::shared_ptr<force> f, _forces) {
    if (((f->part1()->type() == typeAnalyze) or (f->part2()->type() == typeAnalyze)) 
      and (f->part1()->type() != f->part2()->type())) {
      const Eigen::Vector3d radiusVector = rotationAxis.cross(rotationAxis.cross(zeroPoint - f->cP()));
      torqueRet += (radiusVector.cross(f->val())).norm();
    }
  }
  return torqueRet;
};

double snapshot::kinEnergy(int typeAnalyze) {
  double totEnergy = 0.0;
  
  BOOST_FOREACH(std::shared_ptr<particle> p, _particles) {
    if ((p->type() == typeAnalyze or typeAnalyze<0) and (not(p->disabled()))) {
      totEnergy+=p->kinEnergieDouble();
    }
  }
  return totEnergy;
};

double snapshot::potEnergy(int typeAnalyze) {
  double totEnergy = 0.0;
  
   BOOST_FOREACH(std::shared_ptr<force> f, _forces) {
    if ((f->part1()->type() == typeAnalyze) or (f->part1()->type() != f->part2()->type()) or typeAnalyze<0) {
      totEnergy+=f->potEnergyNorm();
    }
  }
  return totEnergy;
};

std::shared_ptr<forceChain> snapshot::forceChainRet() {
  _forceChain = std::make_shared<forceChain>(_particles, _forces);
  return _forceChain;
};

void snapshot::clear() {
  _particles.clear(); _particles.shrink_to_fit();
  _forces.clear(); _forces.shrink_to_fit();
};

snapshot::~snapshot() {
  this->clear();
};
