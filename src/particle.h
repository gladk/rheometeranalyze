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
#include <vector>
#include "interactionori.h"

class particle {
  private:
    Eigen::Vector3d _c, _v, _o;   // Center of mass, linear velocity, angular velocity of the particle
    unsigned long long _id; int _type;     // Particle id, type
    unsigned int _fileId;         // File number of the particle
    double _rad;                  // Particle radius
    double _m;                    // Mass of the particle
    double _d;                    // Density of the particle
    double _volwater;             // Water volume (loaded from data-file)
    
    unsigned short _snapshot;     // Snapshot
    
    Eigen::Matrix3d _axisMatrix;  // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]

    
    bool _disable;                // Disable particle, if it is out of region
    bool _calculateVel;           // Whether the velocity matrix was calculated
    Eigen::Vector3d _vZylindrical;// Linear velocity in zylindrical coordinates (dR, dZ, dF)
    Eigen::Vector3d _posZyl;      // Position in zylindrical coordinates (Rho, Z, phi)
    std::vector <std::shared_ptr<particle> > _contactParticles;     // Vector, where all contacting particles (pointers) are saved
    std::vector <std::shared_ptr<particle> > _contactParticlesWet;  // Vector, where only wet contacting particles (pointers) are saved
    
    std::shared_ptr<interactionori> _normContOri;  // Orientation of normal contacts
    std::shared_ptr<interactionori> _capiContOri;  // Orientation of capillary contacts
    int _sizeIntOri;                               // Size of array of interaction orientations
    
    Eigen::Matrix3d _stressTensor;     // Stress tensor in cylindrical coordinates for normal contacts
    Eigen::Matrix3d _stressTensorCap;  // Stress tensor in cylindrical coordinates for capillar contacts
    double _press, _tau;
    
    bool _shearBand; 
    int _highStressed;
    
  public:
    particle(unsigned long long id, int type, unsigned int fileid, 
             double rad, double mass, double dens, Eigen::Vector3d c, 
             Eigen::Vector3d v, Eigen::Vector3d o, 
             double volwater);
    particle();
    unsigned long long id() {return _id;};
    int type() {return _type;};
    double rad() {return _rad;};
    double mass() {return _m;};
    double density() {return _d;};
    double volwater() {return _volwater;};
    unsigned int  fileId() {return _fileId;};
    double realAngular();    //_v.dot(df)
    Eigen::Vector3d c() {return _c;}
    void c(Eigen::Vector3d cTmp) {_c = cTmp; _calculateVel = false;}
    Eigen::Vector3d v() {return _v;}
    Eigen::Vector3d vZyl();
    Eigen::Vector3d o() {return _o;}
    Eigen::Vector3d dr() {return _axisMatrix.row(0);}
    Eigen::Vector3d dz() {return _axisMatrix.row(1);}
    Eigen::Vector3d df() {return _axisMatrix.row(2);}
    Eigen::Vector3d posZyl() {return _posZyl;}
    void setLocalCoord(Eigen::Vector3d loc, Eigen::Quaternion<double> rotateCC);
    void calculateVel();
    double dist() { return _posZyl(0);};
    double height() { return _posZyl(1);};
    double vol() { return 4.0/3.0*M_PI*_rad*_rad*_rad;};
    double relativeWaterVol();
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
    Eigen::Matrix3d kinEnergie();
    double kinEnergieDouble();
    void addStress(Eigen::Matrix3d addStressTensor);
    void addStressCap(Eigen::Matrix3d addStressTensor);
    void addParticleContact(std::shared_ptr<particle> addParticle);
    void addParticleContactWet(std::shared_ptr<particle> addParticle);
    Eigen::Matrix3d stressTensor();
    Eigen::Matrix3d stressTensorCap();
    Eigen::Matrix3d stressTensorAVG();
    Eigen::Matrix3d stressTensorCapAVG();
    double stressPress();
    double stressTau();
    double stressSigma1();
    double stressSigma3();
    unsigned int contacts();
    unsigned int wetContacts();
    double wetContactsAverageDistance();
    void createIntOri(unsigned short intNumb);
    unsigned short intOriSize() {return _sizeIntOri;}
    InteractionsMatrix normContOri();
    InteractionsMatrix capiContOri();
    bool shearBand() {return _shearBand;}
    void shearBandOn() {_shearBand=true;}
    void shearBandOff() {_shearBand=false;}
    unsigned short snapshot() {return _snapshot;}
    void snapshot(unsigned short snapshot) {_snapshot = snapshot;}
    int highStress() {return _highStressed;}
    void highStress(int highStressed) { _highStressed = highStressed;}
    bool wetContact(std::shared_ptr <particle> p);
    int highStressedContacts();
    bool highStressedContact(std::shared_ptr <particle> p);
};
