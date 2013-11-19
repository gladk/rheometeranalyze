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
  _angles = InteractionsMatrix::Zero(sizeori, sizeori);
  _sizeori = sizeori;
};

void interactionori::clear() {
  _angles.Zero(_sizeori, _sizeori);
};

void interactionori::addinteraction(unsigned short theta, unsigned short psi) {
  if (theta < _sizeori and psi < _sizeori) {
    _angles(theta, psi) += 1;
  } else {
    std::cerr << "Theta = " << theta << "(" << _sizeori << "); Psi = " << psi << " (" << _sizeori << "). Exiting..."<<std::endl; 
    exit (EXIT_FAILURE);
  }
};

void interactionori::addinteraction(Eigen::Vector3d lBranch) {
  lBranch.normalize();
  
  //Rotation around axis df (2)
  Eigen::Quaternion<double> qTheta  = Eigen::Quaternion<double>::FromTwoVectors(Eigen::Vector3d::UnitX(), Eigen::Vector3d(lBranch(0), lBranch(1), 0.0)); qTheta.normalize();
  Eigen::AngleAxis<double>  AaTheta = Eigen::AngleAxis<double>(qTheta);
  double angleAaTheta = AaTheta.angle()*AaTheta.axis()(2);
  
  //Rotation around axis dz (1)
  Eigen::Quaternion<double> qPsi   = Eigen::Quaternion<double>::FromTwoVectors(Eigen::Vector3d::UnitX(), Eigen::Vector3d(lBranch(0), 0.0, lBranch(2))); qPsi.normalize();
  Eigen::AngleAxis<double>  AaPsi = Eigen::AngleAxis<double>(qPsi);
  double angleAaPsi = AaPsi.angle()*AaPsi.axis()(1);
  
  this->addinteraction(floor(2*M_PI/_sizeori*(angleAaTheta + M_PI)), 
                       floor(2*M_PI/_sizeori*(angleAaPsi   + M_PI)));
};

InteractionsMatrix interactionori::interactions() {
  return _angles;
};
