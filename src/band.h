#ifndef BANDCLASS
#define BANDCLASS

#define _USE_MATH_DEFINES

#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "config.h"
#include "particle.h"
#include "force.h"

class band {
  private:
    int _id, _idZ, _idR;                                  // Band ids
    double _dZmin, _dZmax, _dRmin, _dRmax;                // Band minimal and maximal sizes
    long long _partNumb, _forceNumb;                      // Number of particles
    std::vector <boost::shared_ptr<particle> > _allPart;  // Vector of particles;
    std::vector <boost::shared_ptr<force> > _allForces;   // Vector of forces;
    
    double _tau, _tauavg, _vol, _volPart;                 // Results, Tau, Vol, Vol of particles
    double _p, _pavg;                                     // Results, Press
    double _vavg;                                         // Results, angular velocity of particles

  public:
    band(int, int, int, double, double, double, double);
    void addParticle(boost::shared_ptr<particle>);
    void addForce(boost::shared_ptr<force>);
    void calculateValues();
    double TauAVG() {return _tauavg;};
    double PressAVG() {return _pavg;};
    double vol() {return _vol;};
    long long partNumb () {return _partNumb;};
    long long forceNumb () {return _forceNumb;};
    boost::shared_ptr<particle> getPart (long long id) { return _allPart[id];}
    boost::shared_ptr<force> getForc (long long id) { return _allForces[id];}
    int id() {return _id;}
    int idZ() {return _idZ;}
    int idR() {return _idR;}
    double tau() {return _tauavg;}
    double press() {return _pavg;}
    double omega() {return _vavg;}
};

class bandRow {
  private:
    std::vector<boost::shared_ptr<band> > _bandAll;
    boost::shared_ptr<configopt> _cfg; 
    boost::shared_ptr<particleRow> _pRow;
    boost::shared_ptr<forceRow> _fRow;
  public:
    bandRow (boost::shared_ptr<configopt>, boost::shared_ptr<particleRow>, boost::shared_ptr<forceRow>);
    void fillBands();
    int getBandR(double);
    int getBandZ(double);
    void calculateValues();
    boost::shared_ptr<band> getBand(unsigned int id) {return _bandAll[id];}
    unsigned int size() {return _bandAll.size();}
};

#endif
