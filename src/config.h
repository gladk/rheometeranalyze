#ifndef CONFIGCLASS
#define CONFIGCLASS

#include <string>
#include <Eigen/Dense>

class configopt {
  private:
    Eigen::Vector3d _c, _o;   // Center point and origin
    double _Din, _Dout, _H;   // Internal and outer diameters
    int _SecRadial, _SecZ;    // Numer of sections, in which the rheometer will
                              // be divided in radial- and Z-direction 
    Eigen::Vector3d _g;       // Gravity direction

    int _nAt;                 // Number of atom string
    int _nDat;                // Begin of data string

  public:
    configopt(const std::string &);
    Eigen::Vector3d get_c(){return _c;};
    Eigen::Vector3d get_o(){return _o;};
    double Din(){return _Din;};
    double Dout(){return _Dout;};
    double H(){return _H;};
    int SecRadial(){return _SecRadial;};
    int SecZ(){return _SecZ;};
    Eigen::Vector3d get_g(){return _g;};
    int nAt(){return _nAt;};
    int nDat(){return _nDat;};
};

#endif
