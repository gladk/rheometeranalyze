#ifndef BANDCLASS
#define BANDCLASS

#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include "config.h"
#include "particle.h"

class band {
  private:
    int _id, _idZ, _idR, _partNum;           // Band ids, particle number
    double _dZmin, _dZmax, _dRmin, _dRmax;   // Band minimal and maximal sizes

  public:
    band(int, int, int, double, double, double, double);
};

class bandRow {
  private:
    std::vector<boost::shared_ptr<band> > _bandAll;
    boost::shared_ptr<configopt> _cfg; 
    boost::shared_ptr<particleRow> _pRow;
  public:
    bandRow (boost::shared_ptr<configopt>, boost::shared_ptr<particleRow>);
    void fillBands();
};

#endif
