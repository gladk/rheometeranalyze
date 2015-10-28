/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014, 2015 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014, 2015 Anton Gladky <gladky.anton@gmail.com>

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

class bandBase {
  protected:
    bool _cleared = false;
    bool _calculated = false;
    int _id = 0, _idZ = 0, _idR = 0, _idF = 0;            // Band ids
    double _dZmin = 0.,                                   // Band minimal and maximal sizes
           _dZmax = 0.,
           _dRmin = 0.,
           _dRmax = 0.,
           _dFmin = 0.,
           _dFmax = 0.; 
    long long _partNumb = 0;                              // Number of particles
    
    Eigen::Matrix3d _stressTensorAVG = Eigen::Matrix3d::Zero(),             //TODO: Shoud be fixed in accumulator!!!
                    _stressTensorCapAVG = Eigen::Matrix3d::Zero();
    
    double _tau = 0., _tauavg = 0., _vol = 0., _volPart = 0.;  // Results, Tau, Vol, Vol of particles
    double _wetContactsAVG = 0., _wetContactDistanceAVG = 0.;  // Average number of wet contacts and average distance
    double _p = 0., _pavg = 0.;                           // Results, Press (trace), hydrostatic stress
    double  _muAVG = 0.;                                  // Results, mu calculated
    double _vavg = 0., _vavgStDev = 0.;                   // Results, angular velocity of particles, standard deviation
    Eigen::Vector3d _vZylavg = Eigen::Vector3d::Zero();   // Results, average velocity of particles in zylindrical coordinates (dR, dZ, dF)
    double _scherRate = 0.;                               // Results, scherrate
    double _volFraction = 0.;                             // Volume fraction
    double _contactNumAVG = 0.;                           // Contact number per particle, AVG
    double _radAvg = 0.;                                  // Average particle radius
    double _I = 0.;                                       // Schear intensity, dimensionless
    double _densAVG = 0.;                                 // Average particle denisty
    double _eta = 0.;                                     // Viscosity
    double _typeAVG = 0.;                                 // Average "type" of particles
    double _volWaterAVG = 0.;                             // Average water volume of particles
    double _volWaterSUM = 0.;                             // Total water volume of particles
    bool _shearBand=false;                                // True, if the band is in shearband
    std::shared_ptr<configopt> _cfg; 
    
    double _dOmegadR = 0.;                                // dOmega/dr
    double _omega0 = 0.;                                  // omega0, value is needed to get a normalized velocity profile (see Fenistein)
    
    double _gamma=0.0;                                    // local strain, deformations according to Sakaie(5)
    
    InteractionsMatrixD  _normContOri, _capiContOri;      // Interaction orientations
    
    double _d50M=0.0;                                     // Average diameter in respect to mass (volume)

  public:
    bandBase(int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax, std::shared_ptr<configopt> cfg );
    bandBase() {};
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
    double omegaCoefVar();
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
    bool shearBand() { return _shearBand;}
    InteractionsMatrixD normContOri();
    InteractionsMatrixD capiContOri();
    void setdOmegadR(double dOmegadR) {_dOmegadR = dOmegadR;};
    double dOmegadR() const {return _dOmegadR;};
    void omega0(const double omega0) {_omega0 = omega0;};
    double omega0() const {return _omega0;};
    void gamma(const double gamma) {_gamma = gamma;}
    double gamma() const {return _gamma;}
    double d50M() const {return _d50M;}
    const bool cleared() const { return _cleared;}
};
  
class band : public bandBase {
  private:
    std::vector <std::shared_ptr<particle> > _allPart;    // Vector of particles;
  
  public:
    band (int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax, std::shared_ptr<configopt> cfg ) :
      bandBase(id, idZ, idR, idF, dRmin, dRmax, dZmin, dZmax, dFmin, dFmax, cfg ) {};
    
    band (const std::vector<std::shared_ptr<band>> &);
    
    void addParticle(std::shared_ptr<particle>);
    std::shared_ptr<particle> getPart (unsigned long long);
    void calculateValues(int numSnapshots);
    void shearBandOn() { setShearBand(true);}
    void shearBandOff() { setShearBand(false);}
    void setShearBand(const bool shearb);
    void clear();
};
