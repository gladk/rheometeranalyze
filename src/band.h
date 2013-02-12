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
    
    double _tau, _tauavg, _tauLocalAvg, _vol, _volPart;   // Results, Tau, Vol, Vol of particles
    double _p, _pavg, _pLocalAvg;                         // Results, Press
    double _vavg;                                         // Results, angular velocity of particles
    double _volFraction;                                  // Volume fraction
    double _contactNumAVG;                                // Contact number per particle, AVG

  public:
    band(int, int, int, double, double, double, double);
    void addParticle(std::shared_ptr<particle>);
    void addForce(std::shared_ptr<force>);
    void calculateValues();
    double TauAVG() {return _tauavg;};
    double PressAVG() {return _pavg;};
    double vol() {return _vol;};
    double volFraction() {return _volFraction;};
    double contactNumAVG() {return _contactNumAVG;};
    long long partNumb () {return _partNumb;};
    long long forceNumb () {return _forceNumb;};
    std::shared_ptr<particle> getPart (long long id) { return _allPart[id];}
    std::shared_ptr<force> getForc (long long id) { return _allForces[id];}
    int id() {return _id;}
    int idZ() {return _idZ;}
    int idR() {return _idR;}
    double tau() {return _tauavg;}
    double press() {return _pavg;}
    double localPress() {return _pLocalAvg;}
    double omega() {return _vavg;}
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
    unsigned int size() {return _bandAll.size();}
};

#endif
