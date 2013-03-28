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

#include "main.h"
#include "bandRow.h"



bandRow::bandRow (std::shared_ptr<configopt> cfg, std::vector<std::shared_ptr<particleRow>> pRow, std::vector<std::shared_ptr<forceRow>> fRow){
  _cfg =  cfg;
  _pRow = pRow;
  _fRow = fRow;
  int i = 0;
  for (int dz = 0; dz < _cfg->SecZ(); dz++) {
    for (int dr = 0; dr < _cfg->SecRadial(); dr++) {
      std::shared_ptr<band> tmpBand (new band(i, dz, dr, _cfg->Din()/2.0+dr*_cfg->dDr(), _cfg->Din()/2.0+(dr+1)*_cfg->dDr(), _cfg->dDz()*dz, _cfg->dDz()*(dz+1)));
      _bandAll.push_back(tmpBand);
      i++;
    }
  }
  fillBands();
  calculateValues();
};
    
void bandRow::fillBands (){
  //Fill bands with particles
  Eigen::Vector3f O = _cfg->get_c();
  Eigen::Vector3f Z = _cfg->get_o();
  
  //Prepare band-vector
  int i=0;
  for (int z=0; z<_cfg->SecZ(); z++){
    for (int r=0; r<_cfg->SecRadial(); r++){
      double dRmin = _cfg->Din()/2.0 + _cfg->dDr()*r;
      double dRmax = _cfg->Din()/2.0 + _cfg->dDr()*(r+1);
      double dZmin = _cfg->dDz()*z;
      double dZmax = _cfg->dDz()*(z+1);
      std::shared_ptr<band> tmpBand (new band(i, z, r, dRmin, dRmax, dZmin, dZmax));
      _bandAll[i] = tmpBand;
      i++;
    }
  }
  
  //Put particles into band
  long long particleRemoved = 0;
  //Put particles
  for (unsigned int i = 0; i < _pRow.size(); i++)  {
    for (int z = 0; z<_pRow[i]->arraySize(); z++) {
      if (_pRow[i]->particleReal(z)) {
        std::shared_ptr<particle> partTemp = _pRow[i]->getP(z);
        Eigen::Vector3f OP = partTemp->c() - O;     //Vector from center to point
        Eigen::Vector3f OPV = Z.cross(OP);          //Vector, temporal
        OPV.normalize();
        Eigen::Vector3f OPV1 = Z.cross(OPV);        //Vector for projection, Vector Dr
        
        OPV1.normalize();
        double dist = OP.dot(-OPV1);
        partTemp->set_dist(dist);
        double height = OP.dot(Z);
        partTemp->set_height(height);
        partTemp->set_axis(-OPV1, Z, OPV);          //dr, dz, dv
        
        //Define band
        int bR = getBandR(dist);
        int bZ = getBandZ(height);
        
        if (bR>=0 and bZ>=0) {
          int bN = bZ*(_cfg->SecRadial()) + bR;
          _bandAll[bN]->addParticle(partTemp);
        } else {
          _pRow[i]->disable(z);    //Disable and remove particle, if they are out of bands
          particleRemoved ++;
        }
      }
    }
  }
  std::cerr<<particleRemoved<<" particles removed"<<std::endl;
  
 long long forceRemoved = 0;
 //Put forces into band
  for (unsigned int i = 0; i < _fRow.size(); i++)  {
    //Put forces
    for (int z = 0; z<_fRow[i]->arraySize(); z++) {
      std::shared_ptr<force> forceTemp = _fRow[i]->getF(z);
      Eigen::Vector3f OP = forceTemp->cP() - O;   //Vector from center to point
      Eigen::Vector3f OPV = Z.cross(OP);          //Vector, temporal
      OPV.normalize();
      Eigen::Vector3f OPV1 = Z.cross(OPV);        //Vector for projection, Vector Dr
      
      OPV1.normalize();
      double dist = OP.dot(-OPV1);
      forceTemp->set_dist(dist);
      double height = OP.dot(Z);
      forceTemp->set_height(height);
      forceTemp->set_axis(-OPV1, Z, OPV);          //dr, dz, df
      forceTemp->set_dg(_cfg->get_g());
      //Define band
      int bR = getBandR(dist);
      int bZ = getBandZ(height);
      
      if (bR>=0 and bZ>=0) {
        int bN = bZ*(_cfg->SecRadial()) + bR;
        forceTemp->set_band(bR, bZ, bN);
        _bandAll[bN]->addForce(forceTemp);
      } else {
        _fRow[i]->disable(z);    //Disable and remove forces, if they are out of bands
        forceRemoved ++;
      }
    }
  }
  std::cerr<<forceRemoved<<" forces removed"<<std::endl;
  
};


int bandRow::getBandR(double dist) {
  if ((dist>=_cfg->Din()/2.0) and (dist<=_cfg->Dout()/2.0)) {
    return floor((dist - _cfg->Din()/2.0)/_cfg->dDr());
  } else {
    return -1;
  }
};

int bandRow::getBandZ(double height) {
  if ((height>=0) and (height<=_cfg->H())) {
    return floor((height)/_cfg->dDz());
  } else {
    return -1;
  }
};

void bandRow::calculateValues () {
  
  // Common values
  for(unsigned int i=0; i<_bandAll.size(); i++) {
    _bandAll[i]->calculateValues(_cfg->numSnapshot());
  }

  // Scherrate
  for(unsigned int i=1; i<_bandAll.size(); i++) {
    if (_bandAll[i]->idR() > _bandAll[i-1]->idR()) {
      double _shearRateTmp, _shearRateTmpA, _shearRateTmpB;
      //Calculate Scherrate
      
      _shearRateTmpA = (_bandAll[i]->vZyl()(2) - _bandAll[i-1]->vZyl()(2))/
                       (_bandAll[i]->midLinedR() - _bandAll[i-1]->midLinedR() -
                       _bandAll[i-1]->vZyl()(2)/_bandAll[i-1]->midLinedR());
      
      if (i>_cfg->SecRadial()) {
        _shearRateTmpB = (_bandAll[i]->vZyl()(2) - _bandAll[i-_cfg->SecRadial()]->vZyl()(2))/
                         (_bandAll[i]->midLinedZ() - _bandAll[i-_cfg->SecRadial()]->midLinedZ());
      } else {
        _shearRateTmpB = 0.0;
      }
      _shearRateTmp =  0.5*sqrt(_shearRateTmpA*_shearRateTmpA + _shearRateTmpB*_shearRateTmpB);
      _bandAll[i]->set_scherRate(_shearRateTmp);    
      
    }
  }
  
  //Create vector of bandShearZones
  std::shared_ptr<band> maxBandShear;
  std::shared_ptr<band> minBandShear;
  bool maxShearTemp = false;
  bool minShearTemp = false;
  
  for(unsigned int i=1; i<_bandAll.size(); i++) {
    if (_bandAll[i]->idR() > _bandAll[i-1]->idR()) {
      if (_bandAll[i-1]->scherRate() < 0.01 and _bandAll[i]->scherRate() > 0.01 and not(minShearTemp)) {
        minShearTemp = true;
        minBandShear = _bandAll[i-1];
      }
      
      if (_bandAll[i]->scherRate() < 0.01 and _bandAll[i-1]->scherRate() > 0.01 and not(maxShearTemp)) {
        maxShearTemp = true;
        maxBandShear = _bandAll[i];
      }
    } else {
      if (maxShearTemp and minShearTemp) {
        std::shared_ptr<bandShearZone> tmpBandShearZone (new bandShearZone(minBandShear, maxBandShear));
        _bandShearZones.push_back(tmpBandShearZone);
      }
      maxShearTemp = false;
      minShearTemp = false;
    }
  }
  if (maxShearTemp and minShearTemp) {
    std::shared_ptr<bandShearZone> tmpBandShearZone (new bandShearZone(minBandShear, maxBandShear));
    _bandShearZones.push_back(tmpBandShearZone);
  }
};

