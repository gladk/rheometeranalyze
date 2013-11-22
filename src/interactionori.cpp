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

#include "interactionori.h"
#include "math_custom.h"
#include <iostream>

interactionori::interactionori(unsigned short sizeori) {
  _angles = InteractionsMatrix::Zero(sizeori*2, sizeori);
  _sizeori = sizeori;
};

void interactionori::clear() {
  _angles.Zero(_sizeori*2, _sizeori);
};

void interactionori::addinteraction(unsigned short theta, unsigned short psi) {
  if (theta < _sizeori*2 and psi < _sizeori) {
    _angles(theta, psi) += 1;
  } else {
    std::cerr << "Theta = " << theta << "(" << _sizeori*2 << "); Psi = " << psi << " (" << _sizeori << "). Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  }
};

void interactionori::addinteraction(Eigen::Vector3d lBranch) {
  lBranch.normalize();
  
  const Eigen::Vector3d convertA = cart_to_sph(lBranch);
  this->addinteraction(floor(convertA(0)/(M_PI/_sizeori)),    // Theta
                       floor(convertA(1)/(M_PI/_sizeori)));   // Psi
};

InteractionsMatrix interactionori::interactions() {
  return _angles;
};
