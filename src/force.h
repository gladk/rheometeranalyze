#ifndef FORCECLASS
#define FORCECLASS

#include <string>
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <boost/shared_ptr.hpp>

class force {
  private:
    Eigen::Vector3f _pos1, _pos2, _val;   // Pos1, Pos2 and force value
    unsigned long long _pid1, _pid2;      // Particles ids

    double _dist;                 // Distance from axis
    double _height;               // Height on Z-axis
    
    Eigen::Matrix3f _axisMatrix;  // [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]
    
    Eigen::Matrix3f _StressTensor;// [drX, drY, drZ]
                                  // [dzX, dzY, dzZ]
                                  // [dfX, dfY, dfZ]
    
    Eigen::Vector3f _dg;          // Gravity-vector
    int _bandR, _bandZ, _bandN;   // Sections in R-Z directions, section id
    Eigen::Vector3f _cP;          // Contact Point
    
    bool _disable;                // Disable force, if it is out of region
    bool _calculateStressTensor;  // Whether the StressTensor was calculated
    
    
    

  public:
    force(unsigned long long, unsigned long long, Eigen::Vector3f, Eigen::Vector3f, Eigen::Vector3f);
    force();
    int bandR() {return _bandR;};
    int bandZ() {return _bandZ;};
    int bandN() {return _bandN;};
    unsigned long long  pid1() {return _pid1;};
    unsigned long long  pid2() {return _pid2;};
    Eigen::Vector3f pos1() {return _pos1;}
    Eigen::Vector3f pos2() {return _pos2;}
    Eigen::Vector3f val() {return _val;}
    Eigen::Vector3f dr() {return _axisMatrix.row(0);}
    Eigen::Vector3f dz() {return _axisMatrix.row(1);}
    Eigen::Vector3f df() {return _axisMatrix.row(2);}
    Eigen::Vector3f dg() {return _dg;}
    Eigen::Vector3f cP() {return _cP;}
    void set_dist(double dist) {_dist=dist;};
    void set_height(double height) {_height=height;};
    double dist() { return _dist;};
    double height() { return _height;};
    void set_dg(Eigen::Vector3f dg) {_dg=dg;};
    void set_band(int bR, int bZ, int bN) {_bandR=bR; _bandZ=bZ; _bandN=bN;};
    void disable() {_disable=true;};
    void enable() {_disable=false;};
    bool disabled() { return _disable; }
    void set_axis(Eigen::Vector3f dr, Eigen::Vector3f dz, Eigen::Vector3f df);
    double Tau();
    double Press();
    void calculateStressTensor();
};

class forceRow {
  private: 
    std::vector <boost::shared_ptr<force> > _allForce;
    boost::shared_ptr<force> _tmpF;
    long long _realForceNum;
  public:
    forceRow();
    void addF(boost::shared_ptr<force> );
    long long arraySize() {return _allForce.size();};
    long long elementsNum();
    bool forceReal(unsigned long long);
    boost::shared_ptr<force> getF(unsigned long long);
    void disable(unsigned long long);
};
#endif
