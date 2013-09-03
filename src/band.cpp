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
#include <numeric>

band::band(int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax ) {
  _id = id;
  _idZ = idZ;
  _idR = idR;
  _idF = idF;
  _dRmin = dRmin;
  _dRmax = dRmax;
  _dZmin = dZmin;
  _dZmax = dZmax;
  _dFmin = dFmin;
  _dFmax = dFmax;
  _partNumb = 0;
  _tau = 0.0; _tauavg = 0.0;
  _p = 0.0; _pavg = 0.0;
  _vavg = 0.0;
  _vavgStDev = 0.0;
  _vol = M_PI*(dRmax*dRmax - dRmin*dRmin)*(dZmax-dZmin) * ((dFmax-dFmin)/(2*M_PI));
  _volPart = 0.0;
  _volFraction = 0.0;
  _contactNumAVG = 0.0;
  std::vector <std::shared_ptr<particle> > _allPart;
  _stressTensorAVG = _stressTensorAVG.Zero();
  _scherRate = 0.0;
  _muAVG = 0.0;
  _radAvg = 0.0;
  _I = 0.0;
  _densAVG = 0.0;
  _eta = 0.0;
  _vZylavg = _vZylavg.Zero();
  _shearBand = false;
  _wetContactsAVG = 0;
  _wetContactDistanceAVG = 0;
};

void band::addParticle(std::shared_ptr<particle> tmpPart) {
  _allPart.push_back(tmpPart);
  _partNumb ++;
};

void band::set_scherRate(double scherRate) {
  _scherRate = fabs(scherRate);
  /*
  * 
  * The formula (3) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
  * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
  * 
  */ 
  if (_densAVG!=0 and _pavg!=0) {
    _I = _scherRate*(2.0*_radAvg)/(sqrt(_pavg/_densAVG));
  } else {
    _I = .0;
  }
  
  if (_scherRate != 0.0) {
   /*
    * 
    * The formula (1) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
    * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
    * 
    */ 
    _eta = _tauavg/_scherRate;
  }
};

void band::calculateValues (int numSnapshots) {
  _volPart = 0.0;
  _volFraction  = 0.0;
  _contactNumAVG = 0.0;
  _wetContactsAVG = 0.0;
  _wetContactDistanceAVG = 0.0;
  
  Eigen::Matrix3d _totalStressTensor = Eigen::Matrix3d::Zero();
  
  
  unsigned long long i = 0;
  std::vector<double> angVelTmpV;
  std::vector<double> radTMPV;
  std::vector<double> densTMP;
  std::vector<Eigen::Vector3d> velZylTMP;
  
  _tauavg = 0.0;
  _pavg = 0.0;
  for(unsigned long long p=0; p<_allPart.size(); p++) {
    if (not(_allPart[p]->disabled())) {
      angVelTmpV.push_back(_allPart[p]->realAngular());
      velZylTMP.push_back(_allPart[p]->vZyl()*_allPart[p]->vol());
      radTMPV.push_back(_allPart[p]->rad());
      densTMP.push_back(_allPart[p]->density());
      _volPart  += _allPart[p]->vol();
      /*
       * 
       * The formula (15, both parts) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
       * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
       * 
       */
      _totalStressTensor += _allPart[p]->kinEnergie() + _allPart[p]->stressTensor();
      
      _contactNumAVG  += _allPart[p]->contacts();
      _wetContactsAVG += _allPart[p]->wetContacts();
      _wetContactDistanceAVG += _allPart[p]->wetContactsAverageDistance();
      i++;
    }
  }
  
  if (i>0) {
    
     /*
     * 
     * The formula (15) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */
   
    _stressTensorAVG = (_totalStressTensor )/_vol/numSnapshots;
    
    /*
    [S_rr  S_rz  S_rf]
    [S_zr  S_zz  S_zf]
    [S_fr  S_fz  S_ff]
    * 
    * Tau   = sqrt(S_rf*S_rf + S_fz*S_fz)
    * Press = sqrt(S_rr*S_rr + S_zz*S_zz)
    */ 
    
    
    /*
     * 
     * The formula (13) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */ 
    _volFraction  = _volPart/_vol/numSnapshots;
    
    
    _contactNumAVG = _contactNumAVG/i;
    _wetContactsAVG = _wetContactsAVG/i;
    _wetContactDistanceAVG = _wetContactDistanceAVG/i;
    
    _vavg = std::accumulate(angVelTmpV.begin(), angVelTmpV.end(), 0.0) / angVelTmpV.size();
    
    
    /*
     * 
     * The formula (14) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */
    _vZylavg = std::accumulate(velZylTMP.begin(), velZylTMP.end(), Eigen::Vector3d(0,0,0));
    _vZylavg = _vZylavg/_vol/_volFraction/numSnapshots;
    
    double vAVGsq_sum = std::inner_product(angVelTmpV.begin(), angVelTmpV.end(), angVelTmpV.begin(), 0.0);
    _vavgStDev = std::sqrt(vAVGsq_sum / angVelTmpV.size() - _vavg * _vavg);
    
    
    
    
    /*
     * 
     * The formula (--) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */
    _pavg = _stressTensorAVG.trace()/3.0;
    
    
    /*
     * 
     * The formula (17) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */
    _tauavg = sqrt(_stressTensorAVG(2)*_stressTensorAVG(2) + _stressTensorAVG(5)*_stressTensorAVG(5));
    
    _radAvg = std::accumulate(radTMPV.begin(), radTMPV.end(), 0.0) / radTMPV.size();
    _densAVG = std::accumulate(densTMP.begin(), densTMP.end(), 0.0) / densTMP.size();
    
    if (_pavg!= 0.0) {
      _muAVG = _tauavg/_pavg;
    } else {
      _muAVG = 0.0;
    }
  }
};

