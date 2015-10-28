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

#include "band.h"
#include <iostream>
#include <numeric>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>

bandBase::bandBase(int id, int idZ, int idR, int idF, double dRmin, double dRmax, double dZmin, double dZmax, double dFmin, double dFmax, std::shared_ptr<configopt> cfg ) {
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

void bandBase::set_scherRate(double scherRate) {
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
  accumulator_set<double, stats<tag::mean > > acc_volWaterAVG;
  accumulator_set<double, stats<tag::sum  > > acc_volWaterSUM;
  accumulator_set<double, stats<tag::mean > > acc_typeAvg;
  
  
  Eigen::Matrix3d _totalStressTensor = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d _stressTensorCap = Eigen::Matrix3d::Zero();
  
  
  unsigned long long i = 0;
  std::vector<double> angVelTmpV;
  accumulator_set<double, stats<tag::variance > > acc_angVelTmpV;
  accumulator_set<double, stats<tag::mean > > acc_radTMPV;
  accumulator_set<double, stats<tag::mean > > acc_densTMP;
  accumulator_set<double, stats<tag::mean > > acc_d50MTMPV;
  
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
      acc_d50MTMPV(std::pow(_allPart[p]->rad(), 3));
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
      acc_volWaterAVG (_allPart[p]->volwater());
      acc_volWaterSUM (_allPart[p]->volwater());
      
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
    _volWaterAVG = mean(acc_volWaterAVG);
    _volWaterSUM =  sum(acc_volWaterSUM);
    
    _vavg = mean(acc_angVelTmpV);
    _vavgStDev = variance(acc_angVelTmpV);
    
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


    /*
     * Kumar 2014, Acta Mechanica, added after discussion with
     * Prof. Luding 2015/08/11
     *
    _tauavg = sqrt((1./3.)*(std::pow((_stressTensorAVG(0) - _stressTensorAVG(4)),2) +
                            std::pow((_stressTensorAVG(0) - _stressTensorAVG(8)),2) +
                            std::pow((_stressTensorAVG(4) - _stressTensorAVG(8)),2) +
                              6*(std::pow(_stressTensorAVG(1), 2) +
                                  std::pow(_stressTensorAVG(5), 2) +
                                  std::pow(_stressTensorAVG(7), 2))
                            ));
    */
    
    
    /*
     * Implementation from Sudeshna Roy and Thomas Weinhart,
     * similar to ours, 2015/08/12
    
    _tauavg = sqrt(std::pow((_stressTensorAVG(2) + _stressTensorAVG(6)), 2) + std::pow((_stressTensorAVG(5) + _stressTensorAVG(7)), 2))/2;

    */

    _radAvg = mean(acc_radTMPV);
    
    _d50M   = std::pow(mean(acc_d50MTMPV), 1.0/3.0);
    
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

double bandBase::omegaCoefVar() {
if (_vavg) {
    return _vavgStDev/_vavg;
  } else {
    return 0;
  }
}
    
InteractionsMatrixD bandBase::normContOri() {
  return _normContOri;
}

InteractionsMatrixD bandBase::capiContOri() {
  return _capiContOri;
}

double bandBase::omegaNorm() const {
  if (_omega0 != 0) {
    return _vavg/_omega0;
  } else {
    return 0;
  }
}

std::shared_ptr<particle> band::getPart (unsigned long long id) {
  if (id < _allPart.size()) {
    return _allPart[id];
  } else {
    return nullptr;
  }
}

void band::clear() {
  _allPart.clear(); _allPart.shrink_to_fit();
  _cleared = true;
}

band::band(const std::vector<std::shared_ptr<band>> & bV) {
  _id     = bV[0]->_id;
  _idZ    = bV[0]->_idZ;
  _idR    = bV[0]->_idR;
  _idF    = bV[0]->_idF;
  _dZmin  = bV[0]->_dZmin;
  _dZmax  = bV[0]->_dZmax;
  _dRmin  = bV[0]->_dRmin;
  _dRmax  = bV[0]->_dRmax;
  _dFmin  = bV[0]->_dFmin;
  _dFmax  = bV[0]->_dFmax;
  _cfg    = bV[0]->_cfg;
  
  using namespace boost::accumulators;
  typedef accumulator_set< double, features< tag::sum, tag::mean  >, double > ac_d;
  
  ac_d  acc_d; 
  
  for (auto b : bV) { acc_d (b->_tau, weight = b->_partNumb); }; _tau = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_tauavg, weight = b->_partNumb); }; _tauavg = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_vol, weight = b->_partNumb); }; _vol = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_volPart, weight = b->_partNumb); }; _volPart = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_wetContactsAVG, weight = b->_partNumb); }; _wetContactsAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_p, weight = b->_partNumb); }; _p = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_pavg, weight = b->_partNumb); }; _pavg = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_muAVG, weight = b->_partNumb); }; _muAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_vavg, weight = b->_partNumb); }; _vavg = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_vavgStDev, weight = b->_partNumb); }; _vavgStDev = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_scherRate, weight = b->_partNumb); }; _scherRate = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_volFraction, weight = b->_partNumb); }; _volFraction = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_contactNumAVG, weight = b->_partNumb); }; _contactNumAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_radAvg, weight = b->_partNumb); }; _radAvg = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_I, weight = b->_partNumb); }; _I = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_densAVG, weight = b->_partNumb); }; _densAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_eta, weight = b->_partNumb); }; _eta = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_typeAVG, weight = b->_partNumb); }; _typeAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_volWaterAVG, weight = b->_partNumb); }; _volWaterAVG = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_volWaterSUM, weight = b->_partNumb); }; _volWaterSUM = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_dOmegadR, weight = b->_partNumb); }; _dOmegadR = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_omega0, weight = b->_partNumb); }; _omega0 = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_gamma, weight = b->_partNumb); }; _gamma = mean(acc_d); acc_d = ac_d();
  for (auto b : bV) { acc_d (b->_d50M, weight = b->_partNumb); }; _d50M = mean(acc_d); acc_d = ac_d();
  
  // Total number of particles
  long long TotalNumbParticles = 0;
  for (auto b : bV) { TotalNumbParticles+=b->_partNumb;};
  for (auto b : bV) { _vZylavg += b->_vZylavg*(b->_partNumb/TotalNumbParticles);};
  for (auto b : bV) { _stressTensorAVG += b->_stressTensorAVG*(b->_partNumb/TotalNumbParticles);};
  for (auto b : bV) { _stressTensorCapAVG += b->_stressTensorCapAVG*(b->_partNumb/TotalNumbParticles);};
  for (auto b : bV) { _normContOri += b->_normContOri*(b->_partNumb/TotalNumbParticles);};
  for (auto b : bV) { _capiContOri += b->_capiContOri*(b->_partNumb/TotalNumbParticles);};
}
