#ifndef CONFIGCLASS
#define CONFIGCLASS

#include <string>
#include <Eigen/Dense>

class configopt {
  private:
    Eigen::Vector3d _c, _o;   // Center point and origin
    double _Din, _Dout;       // Internal and outer diameters
    int _SecRadial, _SecZ;      // Numer of sections, in which the rheometer will
                              // be divided in radial- and Z-direction 
  public:
    configopt(const std::string &);
    Eigen::Vector3d get_c(){return _c;};
    Eigen::Vector3d get_o(){return _o;};
    double get_Din(){return _Din;};
    double get_Dout(){return _Dout;};
    int get_SecRadial(){return _SecRadial;};
    int get_SecZ(){return _SecZ;};
};

#endif
