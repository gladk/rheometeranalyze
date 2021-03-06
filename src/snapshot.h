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

#pragma once

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp> 

#include "particle.h"
#include "force.h"
#include "forceChain.h"

namespace fs = boost::filesystem;

class snapshot {
  private:
    fs::path _particlesFileName;
    fs::path _forcesFileName;
    unsigned long long _timeStep=0;
    double _time=0;
    unsigned short _id;
    std::vector <std::shared_ptr<particle> > _particles;
    std::vector <std::shared_ptr<force> > _forces;
    std::shared_ptr<forceChain> _forceChain;

  public:
    snapshot(fs::path particlesFileName, fs::path forcesFileName, unsigned long long timeStep);
    void setParticlesFileName(fs::path particlesFileName);
    void setForcesFileName(fs::path forcesFileName);
    void setTimeStep(const unsigned long long timeStep, const double dT);
    unsigned long long timeStep() const;
    unsigned long long size() const;
    double time() const;
    fs::path getParticleFile();
    fs::path getForceFile();
    void addParticle(std::shared_ptr<particle> particleTmp);
    void addForce(std::shared_ptr<force> forceTmp);
    std::vector <std::shared_ptr<particle> > particles();
    std::vector <std::shared_ptr<force> > forces();
    double torque(Eigen::Vector3d rotationAxis, Eigen::Vector3d zeroPoint, int typeAnalyze);
    double kinEnergy(int typeAnalyze);
    double potEnergy(int typeAnalyze);
    unsigned short id() const;
    void id(unsigned short id);
    std::shared_ptr<forceChain> forceChainRet();
    ~snapshot();
    void clear();
};
