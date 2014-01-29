/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014 Anton Gladky <gladky.anton@gmail.com>

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
#include "particle.h"

class force {
  private:
    Eigen::Vector3d _val;   // Force values in cartesian coordinates
    Eigen::Vector3d _cP, _cPZ;    // Contact Point
    
    double _volWater, _distCurr, _distCrit;
    
    std::shared_ptr<particle> _part1;     // Pointer into particle1
    std::shared_ptr<particle> _part2;     // Pointer into particle2

    unsigned int _fileId;           // File number of the particle
     
    Eigen::Vector3d _dg;          // Gravity-vector

    Eigen::Matrix3d _axisMatrix;  // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]

    bool _calculateStressTensor;  // Whether the StressTensor was calculated

  public:
    force(std::shared_ptr<particle> part1, std::shared_ptr<particle> part2, 
          unsigned int fileId, Eigen::Vector3d pos1, Eigen::Vector3d pos2, 
          Eigen::Vector3d val, double volWater, double distCurr, double distCrit);
    force();
    unsigned long long  pid1() {return _part1->id();};
    unsigned long long  pid2() {return _part2->id();};
    std::shared_ptr<particle> part1();
    std::shared_ptr<particle> part2();
    unsigned int  fileId() {return _fileId;};
    Eigen::Vector3d pos1() {return _part1->c();};
    Eigen::Vector3d pos2() {return _part2->c();};
    Eigen::Vector3d pos1Z() {return _part1->posZyl();};
    Eigen::Vector3d pos2Z() {return _part2->posZyl();};
    Eigen::Vector3d val() {return _val;};
    Eigen::Vector3d dg() {return _dg;};
    Eigen::Vector3d cP() {return _cP;};
    Eigen::Vector3d nVec();                                 // Normal vector
    Eigen::Vector3d tVec();                                 // Tangential vector
    double deltaN();                                        // Overlap
    double volWater() {return _volWater;}
    double distCurr() {return _distCurr;}
    double distCrit() {return _distCrit;}
    Eigen::Vector3d forceN();                               // Normal force
    Eigen::Vector3d forceT();                               // Tangential force
    double dist() { return _cPZ(0);};
    double height() { return _cPZ(1);};
    void set_dg(Eigen::Vector3d dg) {_dg=dg;};
    void setLocalCoord(Eigen::Vector3d loc, Eigen::Quaternion<double> rotateCC);
    void calculateStressTensor();
};
