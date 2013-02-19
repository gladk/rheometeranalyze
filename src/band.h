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
    void calculateValues();
    double TauAVG() {return _tauavg;};
    double PressAVG() {return _pavg;};
    double vol() {return _vol;};
    double volFraction() {return _volFraction;};
    double radAvg() {return _radAvg;};
    double contactNumAVG() {return _contactNumAVG;};
    double midLinedR() {return ((_dRmax - _dRmin)/2.0 + _dRmin);};
    double midLinedZ() {return ((_dZmax - _dZmin)/2.0 + _dZmin);};
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

class bandRow {
  private:
    std::vector<std::shared_ptr<band> > _bandAll;
    std::shared_ptr<configopt> _cfg; 
    std::shared_ptr<particleRow> _pRow;
    std::shared_ptr<forceRow> _fRow;
  public:
    bandRow (std::shared_ptr<configopt>, std::shared_ptr<particleRow>, std::shared_ptr<forceRow>);
    void fillBands();
    int getBandR(double);
    int getBandZ(double);
    void calculateValues();
    std::shared_ptr<band> getBand(unsigned int id) {return _bandAll[id];}
    std::shared_ptr<band> getBand(unsigned int idR, unsigned int idZ) {return _bandAll[idZ*_cfg->SecRadial() + idR];}
    unsigned int size() {return _bandAll.size();}
};

#endif
