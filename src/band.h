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

#define _USE_MATH_DEFINES

#include <memory>
#include <Eigen/Dense>
#include "config.h"
#include "particle.h"

class band {
  private:
    int _id, _idZ, _idR, _idF;                            // Band ids
    double _dZmin, _dZmax, _dRmin, _dRmax, _dFmin, _dFmax;// Band minimal and maximal sizes
    long long _partNumb;                                  // Number of particles
    std::vector <std::shared_ptr<particle> > _allPart;    // Vector of particles;
    
    Eigen::Matrix3d _stressTensorAVG, _stressTensorCapAVG;
    
    double _tau, _tauavg, _vol, _volPart;                 // Results, Tau, Vol, Vol of particles
    double _wetContactsAVG, _wetContactDistanceAVG;       // Average number of wet contacts and average distance
    double _p, _pavg;                                     // Results, Press (trace), hydrostatic stress
    double  _muAVG;                                       // Results, mu calculated
    double _vavg, _vavgStDev;                             // Results, angular velocity of particles, standard deviation
    Eigen::Vector3d _vZylavg;                             // Results, average velocity of particles in zylindrical coordinates (dR, dZ, dF)
    double _scherRate;                                    // Results, scherrate
    double _volFraction;                                  // Volume fraction
    double _contactNumAVG;                                // Contact number per particle, AVG
    double _radAvg;                                       // Average particle radius
    double _I;                                            // Schear intensity, dimensionless
    double _densAVG;                                      // Average particle denisty
    double _eta;                                          // Viscosity
    double _typeAVG;                                      // Average "type" of particles
    double _volWaterAVG;                                  // Average water volume of particles
    double _volWaterSUM;                                  // Total water volume of particles
    bool _shearBand;                                      // True, if the band is in shearband
    std::shared_ptr<configopt> _cfg; 
    
    double _dOmegadR=0.0;                                 // dOmega/dr
    double _omega0=0.0;                                   // omega0, value is needed to get a normalized velocity profile (see Fenistein)
    
    double _gamma=0.0;                                    // local strain, deformations according to Sakaie(5)
    
    InteractionsMatrixD  _normContOri, _capiContOri;      // Interaction orientations

  public:
    band(int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax, std::shared_ptr<configopt> cfg );
    void addParticle(std::shared_ptr<particle>);
    void calculateValues(int numSnapshots);
    double TauAVG() {return _tauavg;};
    double PressAVG() {return _pavg;};
    Eigen::Matrix3d TensorAVG() {return _stressTensorAVG;};
    Eigen::Matrix3d TensorCapAVG() {return _stressTensorCapAVG;};
    double vol() {return _vol;};
    double volFraction() {return _volFraction;};
    double radAvg() {return _radAvg;};
    double contactNumAVG() {return _contactNumAVG;};
    double midLinedR() {return ((_dRmax - _dRmin)/2.0 + _dRmin);};
    double midLinedZ() {return ((_dZmax - _dZmin)/2.0 + _dZmin);};
    double midLinedF() {return ((_dFmax - _dFmin)/2.0 + _dFmin);};
    double dR() {return (_dRmax - _dRmin);};
    double dZ() {return (_dZmax - _dZmin);};
    long long partNumb () {return _partNumb;};
    std::shared_ptr<particle> getPart (long long id) { return _allPart[id];}
    int id() const {return _id;}
    int idZ() const {return _idZ;}
    int idR() const {return _idR;}
    int idF() const {return _idF;}
    double tau() {return _tauavg;}
    double press() const {return _pavg;}
    double muAVG() const {return _muAVG;}                            
    double omega() const {return _vavg;}
    double omegaNorm() const;
    double omegaStDev() {return _vavgStDev;}
    void set_scherRate(double );
    void set_I(double I) {_I = I;}
    double scherRate() const { return _scherRate;}
    double I() { return _I;}
    double density() { return _densAVG;}
    double eta() { return _eta;}
    double type() { return _typeAVG;}
    double volwaterAVG() { return _volWaterAVG;}
    double volwaterSUM() { return _volWaterSUM;}
    double wetContactDistanceAVG() { return _wetContactDistanceAVG;}
    double wetContactsAVG() { return _wetContactsAVG;}
    Eigen::Vector3d vZyl() { return _vZylavg;}
    double vDf() { return _vZylavg(2);}
    void shearBandOn() { setShearBand(true);}
    void shearBandOff() { setShearBand(false);}
    bool shearBand() { return _shearBand;}
    InteractionsMatrixD normContOri();
    InteractionsMatrixD capiContOri();
    void setShearBand(const bool shearb);
    void setdOmegadR(double dOmegadR) {_dOmegadR = dOmegadR;};
    double dOmegadR() const {return _dOmegadR;};
    void omega0(const double omega0) {_omega0 = omega0;};
    double omega0() const {return _omega0;};
    void gamma(const double gamma) {_gamma = gamma;}
    double gamma() const {return _gamma;}
};
