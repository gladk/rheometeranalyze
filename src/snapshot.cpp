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

#include "snapshot.h"


snapshot::snapshot(fs::path particlesFileName, fs::path forcesFileName, 
                   unsigned long long timeStep) {
  _particlesFileName = particlesFileName;
  _forcesFileName = forcesFileName;
  _timeStep = timeStep;
};

void snapshot::setParticlesFileName(fs::path particlesFileName) {
  _particlesFileName = particlesFileName;
};

void snapshot::setForcesFileName(fs::path forcesFileName) {
  _forcesFileName = forcesFileName;
};

void snapshot::setTimeStep(unsigned long long timeStep) {
  _timeStep = timeStep;
};

unsigned long long snapshot::timeStep() {
  return _timeStep;
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

std::vector <std::shared_ptr<particle> >  snapshot::particles() {
  return _particles;
};

std::vector <std::shared_ptr<force> >  snapshot::forces() {
  return _forces;
};
