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

#ifndef FORCECLASS
#define FORCECLASS

#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <memory>

class force {
  private:
    Eigen::Vector3f _pos1, _pos2, _val;   // Pos1, Pos2 and force value
    Eigen::Vector3f _valZylindrical;      // Force Values in Zylindrical coordinates (dR, dZ, dF)
    Eigen::Vector3f _cPZylindrical;       // Branch-Vector in Zylindrical coordinates (dR, dZ, dF)
    unsigned long long _pid1, _pid2;      // Particles ids

    double _dist;                 // Distance from axis
    double _height;               // Height on Z-axis
    
    unsigned int _fileId;         // File number of the particle
    
    Eigen::Matrix3f _axisMatrix;  // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]
    
    Eigen::Matrix3f _globalStressTensor;  // [drX, drY, drZ]
                                          // [dzX, dzY, dzZ]
                                          // [dfX, dfY, dfZ]
    
    Eigen::Matrix3f _localStressTensor;
    
    Eigen::Vector3f _dg;          // Gravity-vector
    int _bandR, _bandZ, _bandN;   // Sections in R-Z directions, section id
    Eigen::Vector3f _cP;          // Contact Point
    
    bool _disable;                // Disable force, if it is out of region
    bool _calculateStressTensor;  // Whether the StressTensor was calculated
    
    double _radLen;              // Length of the vector from center to contact point
    
    
    

  public:
    force(unsigned long long, unsigned long long, unsigned  int, Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f);
    force();
    int bandR() {return _bandR;};
    int bandZ() {return _bandZ;};
    int bandN() {return _bandN;};
    unsigned long long  pid1() {return _pid1;};
    unsigned long long  pid2() {return _pid2;};
    unsigned int  fileId() {return _fileId;};
    Eigen::Vector3f pos1() {return _pos1;}
    Eigen::Vector3f pos2() {return _pos2;}
    Eigen::Vector3f val() {return _val;}
    Eigen::Vector3f dr() {return _axisMatrix.row(0);}
    Eigen::Vector3f dz() {return _axisMatrix.row(1);}
    Eigen::Vector3f df() {return _axisMatrix.row(2);}
    Eigen::Vector3f dg() {return _dg;}
    Eigen::Vector3f cP() {return _cP;}
    void set_dist(double dist) {_dist=dist;};
    void set_height(double height) {_height=height;};
    double dist() { return _dist;};
    double height() { return _height;};
    void set_dg(Eigen::Vector3f dg) {_dg=dg;};
    void set_band(int bR, int bZ, int bN) {_bandR=bR; _bandZ=bZ; _bandN=bN;};
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
    void set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df);
    void calculateStressTensor();
    Eigen::Matrix3f localStressTensor();
    Eigen::Matrix3f potEnergie();
};

class forceRow {
  private: 
    std::vector <std::shared_ptr<force> > _allForce;
    std::shared_ptr<force> _tmpF;
    long long _realForceNum;
  public:
    forceRow();
    void addF(std::shared_ptr<force> );
    long long arraySize() {return _allForce.size();};
    long long elementsNum();
    bool forceReal(unsigned long long);
    std::shared_ptr<force> getF(unsigned long long);
    void disable(unsigned long long);
};
#endif
