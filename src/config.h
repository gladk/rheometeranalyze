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
    int _cId;                 // Column of Id
    int _cT;                  // Column of type
    int _cC;                  // Column of center
    int _cV;                  // Column of linear velocity
    int _cO;                  // Column of angular velocity

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
    int cId(){return _cId;};
    int cT(){return _cT;};
    int cC(){return _cC;};
    int cV(){return _cV;};
    int cO(){return _cO;};
};

#endif
