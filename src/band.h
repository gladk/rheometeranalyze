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

#define _USE_MATH_DEFINES

#include <memory>
#include <Eigen/Dense>
#include "config.h"
#include "particle.h"
#include "force.h"

class band {
  private:
    int _id, _idZ, _idR;                                  // Band ids
    double _dZmin, _dZmax, _dRmin, _dRmax;                // Band minimal and maximal sizes
    long long _partNumb, _forceNumb;                      // Number of particles
    std::vector <std::shared_ptr<particle> > _allPart;    // Vector of particles;
    std::vector <std::shared_ptr<force> > _allForces;     // Vector of forces;
    
    Eigen::Matrix3f _stressTensorAVG;
    
    double _tau, _tauavg, _vol, _volPart;                 // Results, Tau, Vol, Vol of particles
    double _p, _pavg;                                     // Results, Press (trace), hydrostatic stress
    double  _muAVG;                                       // Results, mu calculated
    double _vavg, _vavgStDev;                             // Results, angular velocity of particles, standard deviation
    Eigen::Vector3f _vZylavg;                             // Results, average velocity of particles in zylindrical coordinates
    double _scherRate;                                    // Results, scherrate
    double _volFraction;                                  // Volume fraction
    double _contactNumAVG;                                // Contact number per particle, AVG
    double _radAvg;                                       // Average particle radius
    double _I;                                            // Schear intensity, dimensionless
    double _densAVG;                                      // Average particle denisty
    double _eta;                                          // Viscosity
    

  public:
    band(int, int, int, double, double, double, double);
    void addParticle(std::shared_ptr<particle>);
    void addForce(std::shared_ptr<force>);
    void calculateValues(int numSnapshots);
    double TauAVG() {return _tauavg;};
    double PressAVG() {return _pavg;};
    Eigen::Matrix3f TensorAVG() {return _stressTensorAVG;};
    double vol() {return _vol;};
    double volFraction() {return _volFraction;};
    double radAvg() {return _radAvg;};
    double contactNumAVG() {return _contactNumAVG;};
    double midLinedR() {return ((_dRmax - _dRmin)/2.0 + _dRmin);};
    double midLinedZ() {return ((_dZmax - _dZmin)/2.0 + _dZmin);};
    double dR() {return (_dRmax - _dRmin);};
    double dZ() {return (_dZmax - _dZmin);};
    long long partNumb () {return _partNumb;};
    long long forceNumb () {return _forceNumb;};
    std::shared_ptr<particle> getPart (long long id) { return _allPart[id];}
    std::shared_ptr<force> getForc (long long id) { return _allForces[id];}
    int id() {return _id;}
    int idZ() {return _idZ;}
    int idR() {return _idR;}
    double tau() {return _tauavg;}
    double press() {return _pavg;}
    double muAVG() {return _muAVG;}                            
    double omega() {return _vavg;}
    double omegaStDev() {return _vavgStDev;}
    void set_scherRate(double );
    void set_I(double I) {_I = I;}
    double scherRate() { return _scherRate;}
    double I() { return _I;}
    double density() { return _densAVG;}
    double eta() { return _eta;}
    Eigen::Vector3f vZyl() { return _vZylavg;}
    double vDf() { return _vZylavg(2);}
};
