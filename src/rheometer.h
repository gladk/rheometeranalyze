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

#include "band.h"
#include "export.h"
#include "main.h"

namespace fs = boost::filesystem;
using namespace std;
class rheometer {
  private:
    std::shared_ptr<configopt> _cfg;
    
    unsigned long long _particleNum;
    unsigned long long _forceNum;
    std::vector<std::shared_ptr <particleRow>> _particleAll;
    std::vector<std::shared_ptr <forceRow>> _forceRow;
    std::shared_ptr <bandRow>     _bandRow;
    std::shared_ptr <snapshotRow> _snapshots;
  public:
    rheometer(std::shared_ptr<configopt>);
    void loadParticles();
    void loadForces(std::shared_ptr<snapshot> loadSnap);
    void calculateLocalDeformations();
};
