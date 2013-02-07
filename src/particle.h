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
    Eigen::Vector3d _dr;          // Dr-vector
    Eigen::Vector3d _dz;          // Dz-vector
    Eigen::Vector3d _df;          // Df-vector
    int _bandR, _bandZ, _bandN;   // Sections in R-Z directions, section id

  public:
    particle(unsigned long long, int, double, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d);
    particle();
    unsigned long long id() {return _id;};
    int type() {return _type;};
    double rad() {return _rad;};
    int bandR() {return _bandR;};
    int bandZ() {return _bandZ;};
    int bandN() {return _bandN;};
    Eigen::Vector3d c() {return _c;}
    Eigen::Vector3d v() {return _v;}
    Eigen::Vector3d o() {return _o;}
    Eigen::Vector3d dr() {return _dr;}
    Eigen::Vector3d dz() {return _dz;}
    Eigen::Vector3d df() {return _df;}
    void set_dist(double dist) {_dist=dist;};
    double dist() { return _dist;};
    void set_dr(Eigen::Vector3d dr) {_dr=dr;};
    void set_dz(Eigen::Vector3d dz) {_dz=dz;};
    void set_df(Eigen::Vector3d df) {_df=df;};
    void set_band(int bR, int bZ, int bN) {_bandR=bR; _bandZ=bZ; _bandN=bN;};
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
