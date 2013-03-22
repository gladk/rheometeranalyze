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

#ifndef BANDCLASS
#define BANDCLASS

#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <vector>
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
    
    Eigen::Matrix3f _localStressTensorAVG;
    Eigen::Matrix3f _globalStressTensorAVG;
    
    double _tau, _tauavg, _tauLocalAvg, _vol, _volPart;   // Results, Tau, Vol, Vol of particles
    double _p, _pavg, _pLocalAvg;                         // Results, Press (trace), hydrostatic stress
    double _dLocalAvg;                                    // Results, Deviatoric stress, SigmaDev
    double _muLocalAVG, _muGlobAVG;                        // Results, mu calculated
    double _vavg, _vavgStDev;                             // Results, angular velocity of particles, standard deviation
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
    double localPress() {return _pLocalAvg;}                            //Trace, Press
    double dLocalAvg() {return _dLocalAvg;}                             //Deviatoric, SigmaD
    double muLocalAVG() {return _muLocalAVG;}                           //Mu Local
    double muGlobAVG() {return _muGlobAVG;}                             //Mu Global
    double omega() {return _vavg;}
    double omegaStDev() {return _vavgStDev;}
    void set_scherRate(double );
    void set_I(double I) {_I = I;}
    double scherRate() { return _scherRate;}
    double I() { return _I;}
    double density() { return _densAVG;}
    double eta() { return _eta;}
};

class bandShearZone {
  private:
    int _id;
    std::shared_ptr<band> _bandMinPTR;
    std::shared_ptr<band> _bandMaxPTR;
  public:
    bandShearZone(std::shared_ptr<band>, std::shared_ptr<band>);
    double RPOS () {return (_bandMaxPTR->midLinedR() - _bandMinPTR->midLinedR())/2.0 + _bandMinPTR->midLinedR();};
    double ZPOS () {return _bandMaxPTR->midLinedZ();};
    double W () {return (_bandMaxPTR->midLinedR() - _bandMinPTR->midLinedR() - _bandMinPTR->dR());};
};

class bandRow {
  private:
    std::vector<std::shared_ptr<band> > _bandAll;
    std::vector<std::shared_ptr<bandShearZone> > _bandShearZones;
    std::shared_ptr<configopt> _cfg; 
    std::vector<std::shared_ptr<particleRow>> _pRow;
    std::vector<std::shared_ptr<forceRow>> _fRow;
  public:
    bandRow (std::shared_ptr<configopt>, std::vector<std::shared_ptr <particleRow>>, std::vector<std::shared_ptr <forceRow>>);
    void fillBands();
    int getBandR(double);
    int getBandZ(double);
    unsigned int getBandShearZonesSize() {return _bandShearZones.size();};
    std::shared_ptr<bandShearZone> getBandShearZones(int id) {return _bandShearZones[id];};
    void calculateValues();
    std::shared_ptr<band> getBand(unsigned int id) {return _bandAll[id];}
    std::shared_ptr<band> getBand(unsigned int idR, unsigned int idZ) {return _bandAll[idZ*_cfg->SecRadial() + idR];}
    unsigned int size() {return _bandAll.size();}
};

#endif
