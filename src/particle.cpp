#include "particle.h"
#include <iostream>

particle::particle(int id, int type, double rad, Eigen::Vector3d c, Eigen::Vector3d v, Eigen::Vector3d o) {
  _id = id;
  _type = type;
  _rad = rad;
  _c = c;
  _v = v;
  _o = o;
  
};

