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

force::force(std::shared_ptr<particle> part1, std::shared_ptr<particle> part2, unsigned int fileId, Eigen::Vector3f pos1, Eigen::Vector3f pos2, Eigen::Vector3f val) {
  _part1 = part1;
  _part2 = part2;
  _pos1 = pos1;
  _pos2 = pos2;
  _val = val;
  _fileId = fileId;
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
  _calculateStressTensor = false;
  _cP = (pos1-pos2)/2.0 + pos1;
  _stressTensor = _stressTensor.Zero();
  _radLen = ((pos2-pos1)/2.0).norm();
  this->calculateStressTensor();
};

force::force() {
  _fileId = -1;
  _pos1 = Eigen::Vector3f::Zero();
  _pos2 = Eigen::Vector3f::Zero();
  _val = Eigen::Vector3f::Zero();
  _cP = Eigen::Vector3f::Zero();
  _valZylindrical = Eigen::Vector3f::Zero();
  _cPZylindrical = Eigen::Vector3f::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
  _calculateStressTensor = false;
  _stressTensor = _stressTensor.Zero();
  _radLen = 0.0;
};

void force::calculateStressTensor() {
  Eigen::Vector3f l = (_pos1-_pos2)/2.0;
  _stressTensor =  _val*l.transpose();
  _part1->addStress(_stressTensor);
  l = (_pos2-_pos1)/2.0;
  _stressTensor =  (-_val)*l.transpose();
  _part2->addStress(_stressTensor);
  _calculateStressTensor  =  true;
};

Eigen::Matrix3f force::potEnergie() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  return _stressTensor;
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

