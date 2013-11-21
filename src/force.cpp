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

#include "force.h"
#include "math_custom.h"
#include <iostream>

force::force(std::shared_ptr<particle> part1, std::shared_ptr<particle> part2, unsigned int fileId, Eigen::Vector3d pos1, Eigen::Vector3d pos2, Eigen::Vector3d val, double volWater, double distCurr, double distCrit) {
  _part1 = part1;
  _part2 = part2;
  
  _volWater = volWater;
  _distCurr = distCurr;
  _distCrit = distCrit;
  
  if (fabs(pos1.norm()/part1->c().norm() - 1.0) > 0.0001 or fabs(pos2.norm()/part2->c().norm() -1.0 ) > 0.0001 ) {
    std::cerr<<"Particle positions in force and particle files are not the same!"<<std::endl;
    std::cerr<<"pid1 "<<part1->id()<<": ("<<pos1[0]<<" "<<pos1[1]<<" "<<pos1[2]<< ") != (" << part1->c()[0]<<" "<<part1->c()[1]<<" "<<part1->c()[2];
    std::cerr<<"pid2 "<<part2->id()<<": ("<<pos2[0]<<" "<<pos2[1]<<" "<<pos2[2]<< ") != (" << part2->c()[0]<<" "<<part2->c()[1]<<" "<<part2->c()[2]<<");   "<<std::endl;
    
    std::cerr<<"dif1 "<<pos1.norm()/part1->c().norm()<<"; dif2 "<<pos2.norm()/part2->c().norm()<<std::endl<<std::endl;
  }
  
  //Updating particle positions as they are more accurate in fstat
  _part1->c(pos1);
  _part2->c(pos2);
  
  
  _cPZ = Eigen::Vector3d::Zero();
  _val = val;
  _fileId = fileId;
  _calculateStressTensor = false;
  _cP = (pos1-pos2)/2.0 + pos1;
  _axisMatrix = _axisMatrix.Zero();
};

force::force() {
  _fileId = -1;
  _val = Eigen::Vector3d::Zero();
  _cP = Eigen::Vector3d::Zero();
  _cPZ = Eigen::Vector3d::Zero();
  _calculateStressTensor = false;
};

void force::calculateStressTensor() {
  Eigen::Vector3d l;
  Eigen::Matrix3d _stressTensor;
  
  
  
  /*
   * 
   * The formula (15, part2) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
   * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
   * 
   */
  
  
  // Cartesian coordinates

  l = (this->pos1()-this->pos2())/2.0;
  
  // Calculate force into cylindrical coordinates  (projections)
  Eigen::Vector3d FZ = get_cyl_rotated_vector(_val, _axisMatrix);
  
  // Calculate branch-vector into cylindrical coordinates  (projections)
  Eigen::Vector3d lZ = get_cyl_rotated_vector(l, _axisMatrix);
  
  // Calculate the stress tensor
  Eigen::Matrix3d _stressTensorTMP = FZ*lZ.transpose();
  
  _stressTensor = _stressTensorTMP;
  
  
  
  // Add force into the particles
  if (_volWater>0) {
    _part1->addParticleContactWet(_part2);
    _part2->addParticleContactWet(_part1);
    // Set the stress tensor into the particles
    _part1->addStressCap(_stressTensor);
    _part2->addStressCap(_stressTensor);
  } else {
    _part1->addParticleContact(_part2);
    _part2->addParticleContact(_part1);
    // Set the stress tensor into the particles
    _part1->addStress(_stressTensor);
    _part2->addStress(_stressTensor);
  }
  
  _calculateStressTensor  =  true;
};

std::shared_ptr<particle> force::part1() {
  return _part1;
};

std::shared_ptr<particle> force::part2() {
  return _part2;
};

double force::deltaN() {
  return (_part1->rad() + _part2->rad()) - (_part1->c() - _part2->c()).norm();
};

Eigen::Vector3d force::forceN() {
  double valN = this->nVec().dot(_val);
  Eigen::Vector3d forceN = valN*this->nVec();
  return forceN;
};

Eigen::Vector3d force::forceT() {
  Eigen::Vector3d forceT = _val - this->forceN();
  return forceT;
};

Eigen::Vector3d force::nVec() {
  Eigen::Vector3d nVec;
  nVec = _part1->c() - _part2->c();
  nVec.normalize();
  return nVec;
};

Eigen::Vector3d force::tVec() {
  Eigen::Vector3d nReturn = this->forceT();
  nReturn.normalize();
  return nReturn;
};

void force::setLocalCoord(Eigen::Vector3d loc, Eigen::Quaternion<double> rotateCC) {
  _cPZ = cart_to_cyl(rotateCC*loc);
  _axisMatrix = get_axes_coord(_cPZ, rotateCC);
  this->calculateStressTensor();
};

