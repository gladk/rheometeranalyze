#ifndef BANDCLASS
#define BANDCLASS

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
    int _id, _idZ, _idR;                                  // Band ids, particle number
    double _dZmin, _dZmax, _dRmin, _dRmax;                // Band minimal and maximal sizes
    long long _partNumb, _forceNumb;                      // Number of particles
    std::vector <boost::shared_ptr<particle> > _allPart;  // Vector of particles;
    std::vector <boost::shared_ptr<force> > _allForces;   // Vector of forces;

  public:
    band(int, int, int, double, double, double, double);
    void addParticle(boost::shared_ptr<particle>);
    void addForce(boost::shared_ptr<force>);
    void calculateValues();
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
};

#endif
