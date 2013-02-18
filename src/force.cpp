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
  _globalStressTensor = _globalStressTensor.Zero();
  _localStressTensor = _localStressTensor.Zero();
  _radLen = ((pos2-pos1)/2.0).norm();
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
  _globalStressTensor = _globalStressTensor.Zero();
  _localStressTensor = _localStressTensor.Zero();
  _radLen = 0.0;
};

double force::Tau() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  double SigmaF = (_val*_radLen).dot(this->df());
  double SigmaR = (_val*_radLen).dot(this->dr());
  double SigmaT = sqrt(SigmaF*SigmaF + SigmaR*SigmaR);
  
  SigmaT = ((_val*_radLen).cross(this->df())).norm();
  return SigmaT;
};

double force::Press() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  double SigmaP = (_val*_radLen).dot(-this->dz());
  return SigmaP;
};

Eigen::Matrix3f force::localStressTensor() {
  if (not(_calculateStressTensor)) {calculateStressTensor();};
  return _localStressTensor*_radLen;
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
  _globalStressTensor = _axisMatrix.cwiseProduct(forceMatrix);
  
  Eigen::Vector3f lpc; lpc = (_cP - _pos1); lpc.normalize();
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
