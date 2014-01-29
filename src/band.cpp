/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014 Anton Gladky <gladky.anton@gmail.com>

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

#include "band.h"
#include <iostream>
#include <numeric>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>

band::band(int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax, std::shared_ptr<configopt> cfg ) {
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
  _stressTensorCapAVG = _stressTensorCapAVG.Zero();
  _scherRate = 0.0;
  _muAVG = 0.0;
  _radAvg = 0.0;
  _I = 0.0;
  _densAVG = 0.0;
  _eta = 0.0;
  _typeAVG = 0.0;
  _vZylavg = _vZylavg.Zero();
  _shearBand = false;
  _wetContactsAVG = 0;
  _wetContactDistanceAVG = 0;
  _cfg = cfg;
  if (cfg->intOri()>0) {
    _normContOri = InteractionsMatrixD::Zero(cfg->intOri(), cfg->intOri());
    _capiContOri = InteractionsMatrixD::Zero(cfg->intOri(), cfg->intOri());
  } else {
    _normContOri = InteractionsMatrixD::Zero(1, 1);
    _capiContOri = InteractionsMatrixD::Zero(1, 1);
  }
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
    _I = _scherRate*(2.0*_radAvg)/(sqrt(fabs(_pavg/_densAVG)));
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
  
  using namespace boost::accumulators;
  accumulator_set<double, stats<tag::mean > > acc_contactNumAVG;
  accumulator_set<double, stats<tag::mean > > acc_wetContactsAVG;
  accumulator_set<double, stats<tag::mean > > acc_wetContactDistanceAVG;
  accumulator_set<double, stats<tag::mean > > acc_typeAvg;
  
  
  Eigen::Matrix3d _totalStressTensor = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d _stressTensorCap = Eigen::Matrix3d::Zero();
  
  
  unsigned long long i = 0;
  std::vector<double> angVelTmpV;
  accumulator_set<double, stats<tag::mean > > acc_angVelTmpV;
  accumulator_set<double, stats<tag::mean > > acc_radTMPV;
  accumulator_set<double, stats<tag::mean > > acc_densTMP;
  
  
  std::vector<Eigen::Vector3d> velZylTMP;
  
  // Boost::accumulator is not used here, because it returns NULL-value by this typedef
  
  std::vector<InteractionsMatrixD> acc_normContOri;
  std::vector<InteractionsMatrixD> acc_capiContOri;
  
  _tauavg = 0.0;
  _pavg = 0.0;
  for(unsigned long long p=0; p<_allPart.size(); p++) {
    if (not(_allPart[p]->disabled())) {
      
      angVelTmpV.push_back(_allPart[p]->realAngular());
      acc_angVelTmpV(_allPart[p]->realAngular());
      
      velZylTMP.push_back(_allPart[p]->vZyl()*_allPart[p]->vol());
      
      acc_radTMPV(_allPart[p]->rad());
      acc_densTMP(_allPart[p]->density());
      
      _volPart  += _allPart[p]->vol();
      
      /*
       * 
       * The formula (15, both parts) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
       * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
       * 
       */
      _totalStressTensor += _allPart[p]->kinEnergie() + _allPart[p]->stressTensor() + _allPart[p]->stressTensorCap();
      _stressTensorCap +=   _allPart[p]->stressTensorCap();
      
      acc_contactNumAVG(_allPart[p]->contacts());
      acc_wetContactsAVG(_allPart[p]->wetContacts());
      acc_wetContactDistanceAVG (_allPart[p]->wetContactsAverageDistance());
      acc_typeAvg (_allPart[p]->type());
      
      if (_cfg->intOri() > 0) {
        acc_normContOri.push_back(_allPart[p]->normContOri().cast<double>());
        if (_allPart[p]->capiContOri().norm()>0.0) {
          acc_capiContOri.push_back(_allPart[p]->capiContOri().cast<double>());
        }
      }
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
    _stressTensorCapAVG = (_stressTensorCap )/_vol/numSnapshots;
    
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
    
    
    _contactNumAVG = mean(acc_contactNumAVG);
    _wetContactsAVG = mean(acc_wetContactsAVG);
    _wetContactDistanceAVG = mean(acc_wetContactDistanceAVG);
    _typeAVG = mean(acc_typeAvg);
    
    _vavg = mean(acc_angVelTmpV);
    
    /*
     * 
     * The formula (14) in GraMat. Rheology of weakly wetted granular materials - a comparison of experimental and numerical data.
     * Ruediger Schwarze · Anton Gladkyy · Fabian Uhlig · Stefan Luding, 2013
     * 
     */
    
    _vZylavg = std::accumulate(velZylTMP.begin(), velZylTMP.end(), Eigen::Vector3d(0,0,0));
    _vZylavg = _vZylavg/_vol/_volFraction/numSnapshots;
    
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
    
    _radAvg = mean(acc_radTMPV);
    
    _densAVG = mean(acc_densTMP);
    
    if (_cfg->intOri() > 0) {
      if (acc_normContOri.size()>0) {
        _normContOri = std::accumulate(acc_normContOri.begin(), acc_normContOri.end(), _normContOri);
        _normContOri  /= acc_normContOri.size();      // Average contact number in every slot
      } else {
        _normContOri = InteractionsMatrixD::Zero(_cfg->intOri()*2, _cfg->intOri());
      }
      
      if (acc_capiContOri.size()>0) {
        _capiContOri = std::accumulate(acc_capiContOri.begin(), acc_capiContOri.end(), _capiContOri);
        _capiContOri  /= acc_capiContOri.size();      // Average contact number in every slot
      } else {
        _capiContOri = InteractionsMatrixD::Zero(_cfg->intOri()*2, _cfg->intOri());
      }
    }
    
    if (_pavg!= 0.0) {
      _muAVG = _tauavg/_pavg;
    } else {
      _muAVG = 0.0;
    }
  }
};

void band::setShearBand(const bool shearb) {
  for(unsigned long long p=0; p<_allPart.size(); p++) {
    if (not(_allPart[p]->disabled())) {
      if (shearb) {
        _allPart[p]->shearBandOn();
        _shearBand = true;
      } else {
        _allPart[p]->shearBandOff();
        _shearBand = false;
      }
    }
  }
}

InteractionsMatrixD band::normContOri() {
  return _normContOri;
};

InteractionsMatrixD band::capiContOri() {
  return _capiContOri;
};

