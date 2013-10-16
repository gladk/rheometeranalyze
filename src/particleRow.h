/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

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

#include <Eigen/Dense>
#include "particle.h"

class particleRow {
  private: 
    std::vector <std::shared_ptr<particle> > _allPart;
    std::shared_ptr<particle> _tmpP;
    long long _realPartNum;
  public:
    particleRow(long long);
    void addP(std::shared_ptr<particle> );
    long long elementsNum();
    long long arraySize() {return _allPart.size();};
    bool particleReal(long long);
    std::shared_ptr<particle> getP(long long);
    void disable(long long);
    void enable(long long);
};
