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
    unsigned long long _id; int _type;     // Particle id, type
    double _rad;                  // Particle radius

    double _dist;                 // Distance from axis
    double _height;               // Height on Z-axis
    Eigen::Vector3d _dr;          // Dr-vector
    Eigen::Vector3d _dz;          // Dz-vector
    Eigen::Vector3d _df;          // Df-vector
    
    bool _disable;                // Disable particle, if it is out of region

  public:
    particle(unsigned long long, int, double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
    particle();
    unsigned long long id() {return _id;};
    int type() {return _type;};
    double rad() {return _rad;};
    double realAngular() { return _v.dot(_df)/_dist;}
    Eigen::Vector3d c() {return _c;}
    Eigen::Vector3d v() {return _v;}
    Eigen::Vector3d o() {return _o;}
    Eigen::Vector3d dr() {return _dr;}
    Eigen::Vector3d dz() {return _dz;}
    Eigen::Vector3d df() {return _df;}
    void set_dist(double dist) {_dist=dist;};
    void set_height(double height) {_height=height;};
    double dist() { return _dist;};
    double height() { return _height;};
    void set_dr(Eigen::Vector3d dr) {_dr=dr;};
    void set_dz(Eigen::Vector3d dz) {_dz=dz;};
    void set_df(Eigen::Vector3d df) {_df=df;};
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
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
    void disable(long long);
    void enable(long long);
};
#endif
