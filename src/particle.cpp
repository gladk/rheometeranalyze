#include "particle.h"
#include <iostream>

particle::particle(unsigned long long id, int type, double rad, Eigen::Vector3d c, Eigen::Vector3d v, Eigen::Vector3d o) {
  _id = id;
  _type = type;
  _rad = rad;
  _c = c;
  _v = v;
  _o = o;
  
};

particle::particle() {
  _id = -1;
  _type = -1;
  _rad = -1;
  _c = Eigen::Vector3d::Zero();
  _v = Eigen::Vector3d::Zero();
  _o = Eigen::Vector3d::Zero();

};

particleRow::particleRow(long long partN ) {
  boost::shared_ptr<particle> a (new particle());
  _allPart.assign( partN, a ); 
};

void particleRow::addP(boost::shared_ptr<particle> part ) {
  if (_allPart.size()<=part->id()) {
    _allPart.resize(part->id()+1);
  };

 _allPart[part->id()] = part;
};
