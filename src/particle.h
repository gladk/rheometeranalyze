#ifndef PARTICLECLASS
#define PARTICLECLASS

#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>

class particle {
  private:
    Eigen::Vector3d _c, _v, _o;   // Center of mass, linear velocity, angular velocity of the particle
    int _id, _type;               // Particle id, type
    double _rad;                  // Particle radius

    double _dist;                 // Distance from axis

  public:
    particle(unsigned long long, int, double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
    particle();
    unsigned long long id() {return _id;};
    int type() {return _type;};
    double rad() {return _rad;};
    Eigen::Vector3d c() {return _c;}
    Eigen::Vector3d v() {return _v;}
    Eigen::Vector3d o() {return _o;}
    void set_dist(double dist) {_dist=dist;};
    double dist() { return _dist;};
};

class particleRow {
  private: 
    std::vector <boost::shared_ptr<particle> > _allPart;
    boost::shared_ptr<particle> _tmpP;
    long long _realPartNum;
  public:
    particleRow(long long);
    void addP(boost::shared_ptr<particle> );
    long long elementsNum();
    long long arraySize() {return _allPart.size();};
    bool particleReal(long long);
    boost::shared_ptr<particle> getP(long long);
};
#endif
