#include "force.h"

force::force(unsigned long long pid1, unsigned long long pid2, Eigen::Vector3f pos1, Eigen::Vector3f pos2, Eigen::Vector3f val) {
  _pid1 = pid1;
  _pid2 = pid2;
  _pos1 = pos1;
  _pos2 = pos2;
  _val = val;
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
  _calculateStressTensor = false;
  _cP = (pos1-pos2)/2.0 + pos1;
  _axisMatrix = _axisMatrix.Zero();
  _StressTensor = _StressTensor.Zero();
  _localStressTensor = _localStressTensor.Zero();
};

force::force() {
  _pid1 = -1; _pid2 = -1;
  _pos1 = Eigen::Vector3f::Zero();
  _pos2 = Eigen::Vector3f::Zero();
  _val = Eigen::Vector3f::Zero();
  _cP = Eigen::Vector3f::Zero();
  _bandR = -1; _bandZ=-1; _bandN=-1;
  _dist = -1; _height = -1;
  _disable = false;
  _calculateStressTensor = false;
  _axisMatrix = _axisMatrix.Zero();
  _StressTensor = _StressTensor.Zero();
  _localStressTensor = _localStressTensor.Zero();
};

double force::Tau() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  double SigmaR = _val.dot(_axisMatrix.row(2));
  double SigmaZ = _val.dot(_axisMatrix.row(1));
  return sqrt(SigmaR*SigmaR + SigmaZ*SigmaZ);
};

double force::Press() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  double SigmaP = sqrt(_val.dot(_dg)*_val.dot(_dg));
  return SigmaP;
};

Eigen::Matrix3f force::localStressTensor() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  return _localStressTensor;
};

void force::set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df) {
  _axisMatrix = _axisMatrix.Zero();
  dr.normalize(); dz.normalize(); df.normalize();
  _axisMatrix << dr, dz, df;
  _axisMatrix.transposeInPlace();
};

void force::calculateStressTensor() {
  Eigen::Matrix3f forceMatrix; forceMatrix << _val, _val, _val;
  forceMatrix.transposeInPlace();
  _StressTensor = _axisMatrix.cwiseProduct(forceMatrix);
  
  Eigen::Vector3f lpc; lpc = (_pos1-_cP); lpc.normalize();
  _localStressTensor = _val*lpc.transpose();
  
  _calculateStressTensor  =  true;
};

forceRow::forceRow() {
  std::vector <std::shared_ptr<force> > _allForce;
  _realForceNum = 0;
};

void forceRow::addF(std::shared_ptr<force> forc ) {
  _allForce.push_back(forc);
  _realForceNum ++;
};

long long forceRow::elementsNum() {
  return _realForceNum;
};

std::shared_ptr<force> forceRow::getF(unsigned long long id) {
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
