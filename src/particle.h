#ifndef PARTICLECLASS
#define PARTICLECLASS

#include <string>
#include <Eigen/Dense>

class particle {
  private:
    Eigen::Vector3d _c, _v, _o;   // Center of mass, linear velocity, angular velocity of the particle
    int _id;                      // Particle id
    double _rad;                  // Particle radius

  public:
    particle(int, double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
    int id() {return _id;};
    double rad() {return _rad;};
    Eigen::Vector3d c() {return _c;}
    Eigen::Vector3d v() {return _v;}
    Eigen::Vector3d o() {return _o;}
};

#endif
