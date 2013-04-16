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

class particle {
  private:
    Eigen::Vector3f _c, _v, _o;   // Center of mass, linear velocity, angular velocity of the particle
    unsigned long long _id; int _type;     // Particle id, type
    unsigned int _fileId;         // File number of the particle
    double _rad;                  // Particle radius
    double _m;                    // Mass of the particle
    double _d;                    // Density of the particle

    double _dist;                 // Distance from axis
    double _height;               // Height on Z-axis
    
    Eigen::Matrix3f _axisMatrix;  // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]

    Eigen::Matrix3f _velMatrix;   // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]
                                  
    bool _disable;                // Disable particle, if it is out of region
    bool _calculateVel;           // Whether the velocity matrix was calculated
    Eigen::Vector3f _vZylindrical;// Linear velocity in Zylindrical coordinates (dR, dZ, dF)

  public:
    particle(unsigned long long, int, unsigned int, double, double, double, Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f);
    particle();
    unsigned long long id() {return _id;};
    int type() {return _type;};
    double rad() {return _rad;};
    double mass() {return _m;};
    double density() {return _d;};
    unsigned int  fileId() {return _fileId;};
    double realAngular();    //_v.dot(df)
    Eigen::Vector3f c() {return _c;}
    Eigen::Vector3f v() {return _v;}
    Eigen::Vector3f vZyl();
    Eigen::Vector3f o() {return _o;}
    Eigen::Vector3f dr() {return _axisMatrix.row(0);}
    Eigen::Vector3f dz() {return _axisMatrix.row(1);}
    Eigen::Vector3f df() {return _axisMatrix.row(2);}
    void set_dist(double dist) {_dist=dist;};
    void calculateVel();
    void set_height(double height) {_height=height;};
    double dist() { return _dist;};
    double height() { return _height;};
    double vol() { return 4.0/3.0*M_PI*_rad*_rad*_rad;};
    void set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df);
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
    Eigen::Matrix3f kinEnergie();
};
