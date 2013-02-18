#include "particle.h"
#include <iostream>

particle::particle(unsigned long long id, int type, double rad, double mass, double dens, Eigen::Vector3f c, Eigen::Vector3f v, Eigen::Vector3f o) {
  _id = id;
  _type = type;
  _rad = rad;
  _c = c;
  _v = v;
  _o = o;
  _m = mass;
  _d = dens;
  _dist = -1; _height = -1;
  _disable = false;
  _axisMatrix = _axisMatrix.Zero();
  _velMatrix = _velMatrix.Zero();
  _calculateVel = false;
};

particle::particle() {
  _id = -1;
  _type = -1;
  _rad = -1;
  _d = -1;
  _m = -1;
  _c = Eigen::Vector3f::Zero();
  _v = Eigen::Vector3f::Zero();
  _o = Eigen::Vector3f::Zero();
  _dist = -1; _height = -1;
  _disable = false;
  _axisMatrix = _axisMatrix.Zero();
  _velMatrix = _velMatrix.Zero();
  _calculateVel = false;
};

void particle::set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df) {
  _axisMatrix = _axisMatrix.Zero();
  dr.normalize(); dz.normalize(); df.normalize();
  _axisMatrix << dr, dz, df;
  _axisMatrix.transposeInPlace();
  if (not(_calculateVel)) { calculateVel();}
};

void particle::calculateVel() {
  Eigen::Matrix3f velTempMatrix; velTempMatrix << _v, _v, _v;
  velTempMatrix.transposeInPlace();
  _velMatrix = _axisMatrix.cwiseProduct(velTempMatrix);
  _calculateVel = true;
};

double particle::realAngular() {
  if (not(_calculateVel)) { calculateVel();};
  return _v.dot(_axisMatrix.row(2))/_dist;
};

particleRow::particleRow(long long partN ) {
   _tmpP = std::shared_ptr <particle> (new particle());
  _realPartNum = 0;
  _allPart.assign( partN, _tmpP ); 
};

void particleRow::addP(std::shared_ptr<particle> part ) {
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

std::shared_ptr<particle> particleRow::getP(long long id) {
  return _allPart[id];
}

