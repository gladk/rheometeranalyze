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
#include <memory>
#include <vector>

typedef Eigen::Matrix<unsigned short, Eigen::Dynamic, Eigen::Dynamic> InteractionsMatrix;

class interactionori {
  private:
    InteractionsMatrix _angles;     //(Theta, Psi)
    unsigned short _sizeori;
  public:
    interactionori(unsigned short sizeori);
    void clear();
    void addinteraction(unsigned short theta, unsigned short psi); 
    void addinteraction(Eigen::Vector3d lBranch); // Add interaction (force) in zylindrical coordinates (dR, dZ, dF)
    InteractionsMatrix interactions();
};
