#ifndef CONFIGCLASS
#define CONFIGCLASS

#include <string>
#include <Eigen/Dense>

class configopt {
  private:
    Eigen::Vector3f _c, _o;   // Center point and origin
    double _Din, _Dout, _H;   // Internal and outer diameters
    int _SecRadial, _SecZ;    // Numer of sections, in which the rheometer will
                              // be divided in radial- and Z-direction 
    Eigen::Vector3f _g;       // Gravity direction

    int _nAt;                 // Number of atom string
    int _nDat;                // Begin of data string
    int _cId;                 // Column of Id
    int _cT;                  // Column of type
    int _tC;                  // Type of particle, which will be calculated. -1 means calculating of all particles
    int _cC;                  // Column of center
    int _cV;                  // Column of linear velocity
    int _cO;                  // Column of angular velocity
    int _cR;                  // Column of radius
    int _maxC;                // Maximal column number
    void _maxColumnCheck(int, int);    // Check whether the new column is maximal
    void _maxColumnCheckForce(int, int);    // Check whether the new column is maximal for Forces
    
    //Force
    int _fAt;                 // Number of forces string
    int _fDat;                // Begin of force string
    int _cPos1;               // Column of Pos1
    int _cPos2;               // Column of Pos2
    int _cPos1ID;             // Column of particle 1 id
    int _cPos2ID;             // Column of particle 2 id
    int _cForc;               // Column of force value
    int _maxCF;               // Maximal column number for forces
    
  public:
    configopt(const std::string &);
    Eigen::Vector3f get_c(){return _c;};
    Eigen::Vector3f get_o(){return _o;};
    double Din(){return _Din;};
    double Dout(){return _Dout;};
    double H(){return _H;};
    int SecRadial(){return _SecRadial;};
    int SecZ(){return _SecZ;};
    Eigen::Vector3f get_g(){return _g;};
    int nAt(){return _nAt;};
    int nDat(){return _nDat;};
    int cId(){return _cId;};
    int cT(){return _cT;};
    int tC(){return _tC;};
    int cC(){return _cC;};
    int cV(){return _cV;};
    int cO(){return _cO;};
    int cR(){return _cR;};
    int maxC(){return _maxC;};
    int maxCF(){return _maxCF;};
    double dDr(){return (_Dout - _Din)/_SecRadial/2.0;};
    double dDz(){return _H/_SecZ;};
    
    //force
    int fAt(){return _fAt;};
    int fDat(){return _fDat;};
    int cPos1(){return _cPos1;};
    int cPos2(){return _cPos2;};
    int cPos1ID(){return _cPos1ID;};
    int cPos2ID(){return _cPos2ID;};
    int cForc(){return _cForc;};
    
};

#endif
