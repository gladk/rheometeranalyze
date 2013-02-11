#include "force.h"

force::force(unsigned long long pid1, unsigned long long pid2, Eigen::Vector3f pos1, Eigen::Vector3f pos2, Eigen::Vector3f val) {
  _pid1 = pid1;
  _pid2 = pid2;
  _pos1 = pos1;
  _pos2 = pos2;
  _val = val;
  _dr = Eigen::Vector3f::Zero();
  _dz = Eigen::Vector3f::Zero();
  _df = Eigen::Vector3f::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
  _cP = (pos1-pos2)/2.0 + pos1;
};

force::force() {
  _pid1 = -1; _pid2 = -1;
  _pos1 = Eigen::Vector3f::Zero();
  _pos2 = Eigen::Vector3f::Zero();
  _val = Eigen::Vector3f::Zero();
  _dr = Eigen::Vector3f::Zero();
  _dz = Eigen::Vector3f::Zero();
  _df = Eigen::Vector3f::Zero();
  _cP = Eigen::Vector3f::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
};

double force::Tau() {
  double SigmaR = _val.dot(_df);
  double SigmaZ = _val.dot(_dz);
  return sqrt(SigmaR*SigmaR + SigmaZ*SigmaZ);
};

double force::Press() {
  double SigmaP = sqrt(_val.dot(_dg)*_val.dot(_dg));
  return SigmaP;
};

forceRow::forceRow() {
  std::vector <boost::shared_ptr<force> > _allForce;
  _realForceNum = 0;
};

void forceRow::addF(boost::shared_ptr<force> forc ) {
  _allForce.push_back(forc);
  _realForceNum ++;
};

long long forceRow::elementsNum() {
  return _realForceNum;
};

boost::shared_ptr<force> forceRow::getF(unsigned long long id) {
  return _allForce[id];
};

void forceRow::disable(unsigned long long id) {
  if (not(_allForce[id]->disabled())) {
    _allForce[id]->disable();
    _realForceNum--;
  }
};

bool forceRow::forceReal(unsigned long long id) {
  if (_allForce[id]->disabled()) {
    return false;
  } else {
    return true;
  }
};
