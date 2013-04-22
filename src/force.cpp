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

force::force(std::shared_ptr<particle> part1, std::shared_ptr<particle> part2, unsigned int fileId, Eigen::Vector3d pos1, Eigen::Vector3d pos2, Eigen::Vector3d val) {
  _part1 = part1;
  _part2 = part2;
  
  if (pos1 != part1->c() or pos2 != part2->c()) {
    std::cerr<<"Particle positions in force and particle files are not the same!"<<std::endl;
    std::cerr<<"pid1 "<<part1->id()<<": ("<<pos1[0]<<" "<<pos1[1]<<" "<<pos1[2]<< ") != (" << part1->c()[0]<<" "<<part1->c()[1]<<" "<<part1->c()[2]<<std::endl;
    std::cerr<<"pid2 "<<part2->id()<<": ("<<pos2[0]<<" "<<pos2[1]<<" "<<pos2[2]<< ") != (" << part2->c()[0]<<" "<<part2->c()[1]<<" "<<part2->c()[2]<<");   "<<std::endl;
  }

  _valZ = Eigen::Vector3d::Zero();
  _cPZ = Eigen::Vector3d::Zero();
  _val = val;
  _fileId = fileId;
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _disable = false;
  _calculateStressTensor = false;
  _cP = (pos1-pos2)/2.0 + pos1;
  _axisMatrix = _axisMatrix.Zero();
};

force::force() {
  _fileId = -1;
  _val = Eigen::Vector3d::Zero();
  _cP = Eigen::Vector3d::Zero();
  _valZ = Eigen::Vector3d::Zero();
  _cPZ = Eigen::Vector3d::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _disable = false;
  _calculateStressTensor = false;
};

void force::set_axis(Eigen::Vector3d dr, Eigen::Vector3d dz, Eigen::Vector3d df) {
  _axisMatrix = _axisMatrix.Zero();
  dr.normalize(); dz.normalize(); df.normalize();
  _axisMatrix << dr, dz, df;
  _axisMatrix.transposeInPlace();
};

void force::calculateStressTensor() {
  Eigen::Vector3d l;
  Eigen::Matrix3d _stressTensor;
  
  // Cartesian coordinates

  l = (this->pos1()-this->pos2())/2.0;

  Eigen::Matrix3d valTempMatrix; valTempMatrix << _val, _val, _val;
  valTempMatrix.transposeInPlace();
  valTempMatrix = _axisMatrix.cwiseProduct(valTempMatrix);
  
  Eigen::Vector3d FZ = Eigen::Vector3d(valTempMatrix.row(0).sum(),
                                       valTempMatrix.row(1).sum(),
                                       valTempMatrix.row(2).sum());
  
  Eigen::Matrix3d lTempMatrix; lTempMatrix << l, l, l;
  lTempMatrix.transposeInPlace();
  lTempMatrix = _axisMatrix.cwiseProduct(lTempMatrix);
 
  Eigen::Vector3d lZ = Eigen::Vector3d(lTempMatrix.row(0).sum(),
                                       lTempMatrix.row(1).sum(),
                                       lTempMatrix.row(2).sum());

  
  Eigen::Matrix3d _stressTensorTMP = FZ*lZ.transpose();
  
  _stressTensor = _stressTensorTMP;
  
  _part1->addStress(_stressTensor);
  _part2->addStress(_stressTensor);

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

Eigen::Vector3d force::nVec() {
  Eigen::Vector3d nVec;
  nVec = _part1->c() - _part2->c();
  nVec.normalize();
  return nVec;
};

