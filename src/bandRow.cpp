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
#include "math_custom.h"

#ifdef ALGLIB
  #include "interpolation.h"
  void function_cx_1_func(const alglib::real_1d_array &c, const alglib::real_1d_array &x, double &func, void *ptr) 
    {
         func = 0.5 + 0.5*erf((x[0] - c[0])/c[1]);
    }
#endif

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
  Eigen::Vector3d O = _cfg->get_c();
  Eigen::Vector3d Z = _cfg->get_o();
  
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
  
  
  Eigen::Quaternion<double> rotateCCh;   // Rotate coordinate system, hin
  Eigen::Quaternion<double> rotateCCz;   // Rotate coordinate system, zurueck
  
  rotateCCh = rotateCCh.setFromTwoVectors(Eigen::Vector3d(0.0,0.0,1.0), Z);
  
  //Put particles into band
  long long particleRemoved = 0;
  
  for (unsigned int i = 0; i < _pRow.size(); i++)  {
    for (int z = 0; z<_pRow[i]->arraySize(); z++) {
      if (_pRow[i]->particleReal(z)) {
        std::shared_ptr<particle> partTemp = _pRow[i]->getP(z);
        Eigen::Vector3d OP = partTemp->c() - O;     //Vector from center to point
        
        partTemp->setLocalCoord(OP, rotateCCh);

        //Define band
        int bR = getBandR(partTemp->dist());
        int bZ = getBandZ(partTemp->height());
        int bN = bZ*(_cfg->SecRadial()) + bR;
        
        if (bR>=0 and bZ>=0 and bN>=0 and bN<_bandAll.size()) {
          _bandAll[bN]->addParticle(partTemp);
        } else {
          _pRow[i]->disable(z);    // Disable and remove particle, if they are out of bands
          particleRemoved ++;
        }
      }
    }
  }
  std::cout<<particleRemoved<<" particles removed"<<std::endl;
  
  // Calculate stresses
  for (unsigned int i = 0; i < _fRow.size(); i++)  {
    for (int j = 0; j<_fRow[i]->arraySize(); j++) {
      std::shared_ptr<force> forceTemp = _fRow[i]->getF(j);
      
      forceTemp->set_dg(_cfg->get_g());
      Eigen::Vector3d OP = forceTemp->cP() - O;   // Vector from center to contact point
      forceTemp->setLocalCoord(OP, rotateCCh);
    }
  }
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
  for(unsigned int i=1; i<_bandAll.size()-1; i++) {
    if (
       (_bandAll[i+1]->idR() > _bandAll[i-1]->idR()) and 
       (_bandAll[i+1]->idZ() == _bandAll[i-1]->idZ())) {
      double _shearRateTmp, _shearRateTmpA, _shearRateTmpB;
      //Calculate Scherrate
      
      _shearRateTmpA = (_bandAll[i+1]->vZyl()(2) - _bandAll[i-1]->vZyl()(2))/
                       (_bandAll[i+1]->midLinedR() - _bandAll[i-1]->midLinedR() -
                       _bandAll[i]->vZyl()(2)/_bandAll[i]->midLinedR());
      
      if (i>_cfg->SecRadial() and (i+_cfg->SecRadial())<_bandAll.size()) {
        _shearRateTmpB = (_bandAll[i+_cfg->SecRadial()]->vZyl()(2) - _bandAll[i-_cfg->SecRadial()]->vZyl()(2))/
                         (_bandAll[i+_cfg->SecRadial()]->midLinedZ() - _bandAll[i-_cfg->SecRadial()]->midLinedZ());
      } else {
        _shearRateTmpB = 0.0;
      }
      _shearRateTmp =  0.5*sqrt(_shearRateTmpA*_shearRateTmpA + _shearRateTmpB*_shearRateTmpB);
      _bandAll[i]->set_scherRate(_shearRateTmp);    
      
    }
  }
  
  
  //Create vector of Create ShearBands
  
  #ifdef ALGLIB
    
    double omega0 = this->getBand(_cfg->SecRadial()-1,0)->omega(); 
    for(int h=0; h<_cfg->SecZ(); h++) {
      alglib::real_2d_array x;
      alglib::real_1d_array y;
      alglib::real_1d_array c = "[0.08, 0.0075]";
      x.setlength(_cfg->SecRadial(), 1);
      y.setlength(_cfg->SecRadial());
      
      for(int r=0; r<_cfg->SecRadial(); r++) {
      
        double omega = this->getBand(r,h)->omega(); 
        double valT = omega/omega0;
        x(r,0) = (this->getBand(r,0))->midLinedR();
        y(r) = valT;
      }
      double epsf = 0;
      double epsx = 0.0000001;
      alglib::ae_int_t maxits = 0;
      alglib::ae_int_t info;
      alglib::lsfitstate state;
      alglib::lsfitreport rep;
      double diffstep = 0.000001;
      lsfitcreatef(x, y, c, diffstep, state);
      lsfitsetcond(state, epsf, epsx, maxits);
      alglib::lsfitfit(state, function_cx_1_func);
      lsfitresults(state, info, c, rep);
      
      const double Rz = c(0);
      const double W = c(1);
      _shearBands.push_back(Eigen::Vector3d(Rz, W, this->getBand(0,h)->midLinedZ()));
      for(int r=0; r<_cfg->SecRadial(); r++) {
        double const Rt = this->getBand(r,h)->midLinedR();
        if ((Rt>=(Rz-W)) and (Rt<=(Rz+W))) {
          this->getBand(r,h)->shearBandOn();
        }
      }
    }
  #endif
};

