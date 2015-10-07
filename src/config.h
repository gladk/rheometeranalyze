/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014, 2015 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014, 2015 Anton Gladky <gladky.anton@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    RheometerAnalyze is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with RheometerAnalyze.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#define BOOST_FILESYSTEM_NO_DEPRECATED
#include <boost/filesystem.hpp> 
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <iostream>
#include <string>
#include <Eigen/Dense>

#include "snapshotRow.h"

class configopt {
  private:
    Eigen::Vector3d _c, _o;   // Center point and origin
    double _Din, _Dout, _H;   // Internal and outer diameters
    int _SecRadial, _SecZ, 
                    _SecFi;   // Number of sections, in which the rheometer will
                              // be divided in radial-, Z- and Fi-directions
    Eigen::Vector3d _g;       // Gravity direction
      
    double _dT;               // Delta time
    int _nAt;                 // Number of atom string
    int _nPSt;                // Number of timestep string (particles)
    int _nDat;                // Begin of data string
    int _cId;                 // Column of Id
    int _cT;                  // Column of type
    int _tC;                  // Type of particle, which will be calculated. -1 means calculating of all particles
    int _tR;                  // Type of particles, which rotates. Needs for calculation of torque
    int _cC;                  // Column of center
    int _cV;                  // Column of linear velocity
    int _cO;                  // Column of angular velocity
    int _cR;                  // Column of radius
    int _cM;                  // Column of mass
    int _cD;                  // Column of density
    int _cVolWaterP;          // Column of volWater value (for Particle)
    int _maxC;                // Maximal column number
    void _maxColumnCheck(int, int);    // Check whether the new column is maximal
    void _maxColumnCheckForce(int, int);    // Check whether the new column is maximal for Forces
    std::string _FOutput;          //Folder for output files
    double _aS, _aE;           // Angle to analyze in degrees, aStart and aEnd
    double _aDiff;             // Range in degrees
    
    // Forces
    int _fAt;                 // Number of forces string
    int _nFSt;                // Number of timestep string (forces)
    int _fDat;                // Begin of force string
    int _tF;                  // Type of particles, which will be calculated (Forces). -1 means calculating of all particles
    int _cPos1;               // Column of Pos1
    int _cPos2;               // Column of Pos2
    int _cPos1ID;             // Column of particle 1 id
    int _cPos2ID;             // Column of particle 2 id
    int _cForc;               // Column of force value
    int _cVolWater;           // Column of volWater value
    int _cDistCurr;           // Column of cDistCurr value
    int _cDistCrit;           // Column of cDistCrit value
    int _maxCF;               // Maximal column number for forces
    int _wetParticle;         // Number of slots for the particle wet analyze
    
    int _intOri;              // Number of slots for analyze the orientation of interactions
  
    // Others
    std::shared_ptr<snapshotRow> _snapshotRow;    // Row of snapshots
    unsigned short  _vtk;      // True, if VTK-file will be created (0 - no export, 1 - all data to export, 2 - only band-data to export)
    bool  _utwente;            // True, if UTwente-files will be created
    bool  _contact;            // True, if contact-analyze should be performed
    bool  _followContact;      // True, if contact-follow-analyze should be performed
    bool  _noShearBand=false;  // True, if we do not need shear bands
    bool  _increment=false;    // True, if incremental calculation is needed
    
    double _omega0;            // Rotating velocity
    
  public:
    configopt(const std::string &);
    Eigen::Vector3d get_c(){return _c;};
    Eigen::Vector3d get_o(){return _o;};
    double Din() const {return _Din;};
    double Dout() const {return _Dout;};
    double H()  const {return _H;};
    double aS() const {return _aS;};
    double aE() const {return _aE;};
    double aDiff() const {return _aDiff;};
    unsigned int SecRadial(){return _SecRadial;};
    unsigned int SecZ(){return _SecZ;};
    unsigned int SecFi(){return _SecFi;};
    Eigen::Vector3d get_g(){return _g;};
    int nAt() const {return _nAt;};
    int nPSt() const {return _nPSt;};
    int nDat() const {return _nDat;};
    int cId() const {return _cId;};
    int cT() const {return _cT;};
    int tC() const {return _tC;};
    int tR() const {return _tR;};
    int cC() const {return _cC;};
    int cV() const {return _cV;};
    int cO() const {return _cO;};
    int cR() const {return _cR;};
    int cD() const {return _cD;};
    int cM() const {return _cM;};
    int cVolWaterP() const {return _cVolWaterP;};
    int maxC() const {return _maxC;};
    int maxCF() const {return _maxCF;};
    int numSnapshot() const {return _snapshotRow->size();};
    void setSnapshot(std::shared_ptr<snapshotRow> snapshotRowTMP) {_snapshotRow = snapshotRowTMP;}
    void unSetSnapshot() {_snapshotRow = nullptr;}
    std::shared_ptr<snapshotRow> snapshot() {return _snapshotRow;}
    double dDr(){return (_Dout - _Din)/_SecRadial/2.0;};
    double dDz(){return _H/_SecZ;};
    double dDf(){return (_aE-_aS)/_SecFi;};
    double dT(){return _dT;};
    std::string FOutput(){return _FOutput;};
    void FOutput(std::string FOutput){_FOutput=FOutput;};
    
    // Force
    int fAt(){return _fAt;};
    int nFSt(){return _nFSt;};
    int tF(){return _tF;};
    int fDat(){return _fDat;};
    int cPos1(){return _cPos1;};
    int cPos2(){return _cPos2;};
    int cPos1ID(){return _cPos1ID;};
    int cPos2ID(){return _cPos2ID;};
    int cForc(){return _cForc;};
    int cVolWater(){return _cVolWater;};
    int cDistCurr(){return _cDistCurr;};
    int cDistCrit(){return _cDistCrit;};
    
    // Others
    void setVtk(unsigned short setVtk) {_vtk=setVtk;}
    void unSetVtk() {_vtk=0;}
    unsigned short Vtk() {return _vtk;}
    
    void setUtwente() {_utwente=true;}
    void unSetUtwente() {_utwente=false;}
    bool Utwente() {return _utwente;}  

    void setContact() {_contact=true;}
    void unSetContact() {_contact=false;}
    bool contact() {return _contact;} 
    
    void setFollowContact() {_followContact=true;} 
    void unSetFollowContact() {_followContact=false;} 
    bool followContact() {return _followContact;}  
    
    void setWetParticle(int i = 0) {_wetParticle=i;} 
    void unSetWetParticle() {_wetParticle=-1;} 
    int wetParticle() {return _wetParticle;}  
    
    void setIntOri(int intOri) {_intOri=intOri;}
    void unSetIntOri() {_intOri=-1;}
    int intOri() {return _intOri;}
    
    void setOmega0(double omega0) {_omega0=omega0;}
    void unSetOmega0() {_omega0=0;}
    double omega0() const {return _omega0;}
    
    void setNoShearBand() {_noShearBand=true;}
    void unSetNoShearBand() {_noShearBand=false;}
    bool noShearBand() const {return _noShearBand;}
    
    void setIncrement() {_increment=true;}
    void unSetIncrement() {_increment=false;}
    bool increment() const {return _increment;}
    
    // vector, containing local deformations, gamma
    std::vector<double> _gamma;
    // current time step
    double _timeCur = -1.0;
    unsigned long long _timeStepCur = 0;
};
