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

#include "force.h"
#include <iostream>

force::force(std::shared_ptr<particle> part1, std::shared_ptr<particle> part2, unsigned int fileId, Eigen::Vector3f pos1, Eigen::Vector3f pos2, Eigen::Vector3f val) {
  _part1 = part1;
  _part2 = part2;
  
  if (pos1 != part1->c() or pos2 != part2->c()) {
    std::cerr<<"Particle positions in force and particle files are not the same1!"<<std::endl;
    std::cerr<<part1->id()<< "pid1, "<<part2->id()<< "pid2"<<std::endl<<std::endl;
  }

  _valZ = Eigen::Vector3f::Zero();
  _cPZ = Eigen::Vector3f::Zero();
  _val = val;
  _fileId = fileId;
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _disable = false;
  _calculateStressTensor = false;
  _cP = (pos1-pos2)/2.0 + pos1;
  this->calculateStressTensor();
};

force::force() {
  _fileId = -1;
  _val = Eigen::Vector3f::Zero();
  _cP = Eigen::Vector3f::Zero();
  _valZ = Eigen::Vector3f::Zero();
  _cPZ = Eigen::Vector3f::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _disable = false;
  _calculateStressTensor = false;
};

void force::calculateStressTensor() {
  Eigen::Vector3f l;
  Eigen::Matrix3f _stressTensor;
  
  // Cartesian coordinates
  
  l = (this->pos1()-this->pos2())/2.0;
  _stressTensor =  _val*l.transpose();
  _part1->addStress(_stressTensor);
  
  l = (this->pos2()-this->pos1())/2.0;
  _stressTensor =  (-_val)*l.transpose();
  _part2->addStress(_stressTensor);
  
  
  // Cylindrical coordinates
  /*
  l = (this->pos1Z()-this->pos2Z())/2.0;
  _stressTensor =  _valZ*l.transpose();
  _part1->addStress(_stressTensor);
  
  l = (this->pos2Z()-this->pos1Z())/2.0;
  _stressTensor =  (-_valZ)*l.transpose();
  _part2->addStress(_stressTensor);
  */
  _calculateStressTensor  =  true;
};

std::shared_ptr<particle> force::part1() {
  return _part1;
};

std::shared_ptr<particle> force::part2() {
  return _part2;
};

double force::deltaN() {
  return (_part1->c() - _part2->c()).norm() - _part1->rad() + _part2->rad();
};

double force::valN() {
  double valN = this->nVec().dot(_val);
  return fabs(valN);
};

Eigen::Vector3f force::nVec() {
  Eigen::Vector3f nVec;
  nVec = _part1->c() - _part2->c();
  nVec.normalize();
  return nVec;
};

