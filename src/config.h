/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Foobar is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONFIGCLASS
#define CONFIGCLASS

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp> 
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <iostream>
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
    int _cM;                  // Column of mass
    int _cD;                  // Column of density
    int _maxC;                // Maximal column number
    void _maxColumnCheck(int, int);    // Check whether the new column is maximal
    void _maxColumnCheckForce(int, int);    // Check whether the new column is maximal for Forces
    std::string _FOutput;          //Folder for output files
    
    //Force
    int _fAt;                 // Number of forces string
    int _fDat;                // Begin of force string
    int _cPos1;               // Column of Pos1
    int _cPos2;               // Column of Pos2
    int _cPos1ID;             // Column of particle 1 id
    int _cPos2ID;             // Column of particle 2 id
    int _cForc;               // Column of force value
    int _maxCF;               // Maximal column number for forces
    
    //other
    int _numSnapshot;          // Number of snapshots
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
    int cD(){return _cD;};
    int cM(){return _cM;};
    int maxC(){return _maxC;};
    int maxCF(){return _maxCF;};
    int numSnapshot(){return _numSnapshot;};
    void setSnapshot(int numSnapshot) {_numSnapshot = numSnapshot;}
    double dDr(){return (_Dout - _Din)/_SecRadial/2.0;};
    double dDz(){return _H/_SecZ;};
    std::string FOutput(){return _FOutput;};
    
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
