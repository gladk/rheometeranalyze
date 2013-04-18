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

#pragma once

#include <Eigen/Dense>
#include <memory>
#include "particle.h"

class force {
  private:
    Eigen::Vector3f _val, _valZ;   // Force values in cartesian and Zylindrical coordinates (dR, dZ, dF)
    Eigen::Vector3f _cP,  _cPZ;    // Contact Point
    
    std::shared_ptr<particle> _part1;     // Pointer into particle1
    std::shared_ptr<particle> _part2;     // Pointer into particle2

    unsigned int _fileId;           // File number of the particle
     
    Eigen::Vector3f _dg;          // Gravity-vector
    int _bandR, _bandZ, _bandN;   // Sections in R-Z directions, section id
    
    
    bool _disable;                // Disable force, if it is out of region
    bool _calculateStressTensor;  // Whether the StressTensor was calculated
    

  public:
    force(std::shared_ptr<particle>, std::shared_ptr<particle>, unsigned  int, Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f);
    force();
    int bandR() {return _bandR;};
    int bandZ() {return _bandZ;};
    int bandN() {return _bandN;};
    unsigned long long  pid1() {return _part1->id();};
    unsigned long long  pid2() {return _part2->id();};
    std::shared_ptr<particle> part1();
    std::shared_ptr<particle> part2();
    unsigned int  fileId() {return _fileId;};
    Eigen::Vector3f pos1() {return _part1->c();};
    Eigen::Vector3f pos2() {return _part2->c();};
    Eigen::Vector3f pos1Z() {return _part1->posZyl();};
    Eigen::Vector3f pos2Z() {return _part2->posZyl();};
    Eigen::Vector3f val() {return _val;};
    Eigen::Vector3f dg() {return _dg;};
    Eigen::Vector3f cP() {return _cP;};
    Eigen::Vector3f nVec();
    double deltaN();                                        // Overlap
    double valN();                                          // Normal force
    double dist() { return _cPZ(0);};
    double height() { return _cPZ(1);};
    void set_dg(Eigen::Vector3f dg) {_dg=dg;};
    void set_band(int bR, int bZ, int bN) {_bandR=bR; _bandZ=bZ; _bandN=bN;};
    void set_cPZ(Eigen::Vector3f cPZ) {_cPZ = cPZ;};
    void set_valZ(Eigen::Vector3f valZ) {_valZ = valZ;};
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
    void calculateStressTensor();
};
