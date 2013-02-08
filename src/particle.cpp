#include "particle.h"
#include <iostream>

particle::particle(unsigned long long id, int type, double rad, Eigen::Vector3d c, Eigen::Vector3d v, Eigen::Vector3d o) {
  _id = id;
  _type = type;
  _rad = rad;
  _c = c;
  _v = v;
  _o = o;
  _dr = Eigen::Vector3d::Zero();
  _dz = Eigen::Vector3d::Zero();
  _df = Eigen::Vector3d::Zero();
  _dist = -1; _height = -1;
  _disable = false;
};

particle::particle() {
  _id = -1;
  _type = -1;
  _rad = -1;
  _c = Eigen::Vector3d::Zero();
  _v = Eigen::Vector3d::Zero();
  _o = Eigen::Vector3d::Zero();
  _dr = Eigen::Vector3d::Zero();
  _dz = Eigen::Vector3d::Zero();
  _df = Eigen::Vector3d::Zero();
  _dist = -1; _height = -1;
  _disable = false;
};

particleRow::particleRow(long long partN ) {
   _tmpP = boost::shared_ptr <particle> (new particle());
  _realPartNum = 0;
  _allPart.assign( partN, _tmpP ); 
};

void particleRow::addP(boost::shared_ptr<particle> part ) {
  if (_allPart.size()<=part->id()) {
    _allPart.resize(part->id()+1, _tmpP);
  };
  if (_allPart[part->id()] == _tmpP) {
    _allPart[part->id()] = part;
    _realPartNum ++;
  } else {
    std::cerr<<"The particles have equal IDs. "<<part->id()<<" Aborting."<<std::endl;
    exit (EXIT_FAILURE);
  }
};

long long particleRow::elementsNum() {
  return _realPartNum;
};

bool particleRow::particleReal(long long id) {
  if (_allPart[id] == _tmpP or _allPart[id]->disabled()) {
    return false;
  } else {
    return true;
  }
};

void particleRow::disable(long long id) {
  if (not(_allPart[id]->disabled())) {
    _allPart[id]->disable();
    _realPartNum--;
  }
}

boost::shared_ptr<particle> particleRow::getP(long long id) {
  return _allPart[id];
}

