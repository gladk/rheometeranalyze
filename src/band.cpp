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

#include "band.h"
#include <iostream>

band::band(int id, int idZ, int idR, double dRmin, double dRmax, double dZmin, double dZmax) {
  _id = id;
  _idZ = idZ;
  _idR = idR;
  _dRmin = dRmin;
  _dRmax = dRmax;
  _dZmin = dZmin;
  _dZmax = dZmax;
  _partNumb = 0;
  _forceNumb = 0;
  _tau = 0.0; _tauavg = 0.0;
  _p = 0.0; _pavg = 0.0;
  _vavg = 0.0;
  _vavgStDev = 0.0;
  _vol = M_PI*(dRmax*dRmax - dRmin*dRmin)*(dZmax-dZmin);
  _volPart = 0.0;
  _volFraction = 0.0;
  _contactNumAVG = 0.0;
  std::vector <std::shared_ptr<particle> > _allPart;
  std::vector <std::shared_ptr<force> > _allForces;
  _stressTensorAVG = _stressTensorAVG.Zero();
  _scherRate = 0.0;
  _muAVG = 0.0;
  _radAvg = 0.0;
  _I = 0.0;
  _densAVG = 0.0;
  _eta = 0.0;
  _vZylavg = _vZylavg.Zero();
};

void band::addParticle(std::shared_ptr<particle> tmpPart) {
  _allPart.push_back(tmpPart);
  _partNumb ++;
};

void band::addForce(std::shared_ptr<force> tmpForc) {
  _allForces.push_back(tmpForc);
  _forceNumb ++;
};

void band::set_scherRate(double scherRate) {
  _scherRate = fabs(scherRate);
  _I = _scherRate*(2.0*_radAvg)/(sqrt(_densAVG/_pavg));
  _eta = _tauavg/_scherRate;
};

void band::calculateValues (int numSnapshots) {
  _volPart = 0.0;
  _volFraction  = 0.0;
  
  Eigen::Matrix3f _totalStressTensor = Eigen::Matrix3f::Zero();
  
  
  unsigned long long i = 0;
  std::vector<double> angVelTmpV;
  std::vector<double> radTMPV;
  std::vector<double> densTMP;
  std::vector<Eigen::Vector3f> velZylTMP;
  
  _tauavg = 0.0;
  _pavg = 0.0;
  for(unsigned long long p=0; p<_allPart.size(); p++) {
    if (not(_allPart[p]->disabled())) {
      angVelTmpV.push_back(_allPart[p]->realAngular());
      velZylTMP.push_back(_allPart[p]->vZyl());
      radTMPV.push_back(_allPart[p]->rad());
      densTMP.push_back(_allPart[p]->density());
      _volPart  += _allPart[p]->vol();
      _totalStressTensor += _allPart[p]->kinEnergie() + _allPart[p]->stressTensor();
      i++;
    }
  }
  
  if (i>0) {
    
    _stressTensorAVG = (_totalStressTensor )/_vol/numSnapshots;
    
    /*
    [S_rr  S_rz  S_rf]
    [S_zr  S_zz  S_zf]
    [S_fr  S_fz  S_ff]
    * 
    * Tau   = sqrt(S_rf*S_rf + S_fz*S_fz)
    * Press = sqrt(S_rr*S_rr + S_zz*S_zz)
    */ 
    
    _volFraction  = _volPart/_vol/numSnapshots;
    _contactNumAVG = (double)_allForces.size()/i;
    _vavg = std::accumulate(angVelTmpV.begin(), angVelTmpV.end(), 0.0) / angVelTmpV.size();
    
    _vZylavg = std::accumulate(velZylTMP.begin(), velZylTMP.end(), Eigen::Vector3f(0,0,0)) / velZylTMP.size();
    
    double vAVGsq_sum = std::inner_product(angVelTmpV.begin(), angVelTmpV.end(), angVelTmpV.begin(), 0.0);
    _vavgStDev = std::sqrt(vAVGsq_sum / angVelTmpV.size() - _vavg * _vavg);
    
    _pavg = sqrt(_stressTensorAVG.row(1).norm()*_stressTensorAVG.row(1).norm() + _stressTensorAVG.row(0).norm()*_stressTensorAVG.row(0).norm());
    _tauavg = _stressTensorAVG.row(2).norm();;
    
    _radAvg = std::accumulate(radTMPV.begin(), radTMPV.end(), 0.0) / radTMPV.size();
    _densAVG = std::accumulate(densTMP.begin(), densTMP.end(), 0.0) / densTMP.size();
    
    if (_pavg!= 0.0) {
      _muAVG = _tauavg/_pavg;
    } else {
      _muAVG = 0.0;
    }
  }
};

