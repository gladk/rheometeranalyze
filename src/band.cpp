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
  _vol = M_PI/4.0*(dRmax*dRmax - dRmin*dRmin)*(dZmax-dZmin);
  std::vector <boost::shared_ptr<particle> > _allPart;
  std::vector <boost::shared_ptr<force> > _allForces;
};

void band::addParticle(boost::shared_ptr<particle> tmpPart) {
  _allPart.push_back(tmpPart);
  _partNumb ++;
};

void band::addForce(boost::shared_ptr<force> tmpForc) {
  _allForces.push_back(tmpForc);
  _forceNumb ++;
};


bandRow::bandRow (boost::shared_ptr<configopt> cfg, boost::shared_ptr<particleRow> pRow, boost::shared_ptr<forceRow> fRow){
  _cfg =  cfg;
  _pRow = pRow;
  _fRow = fRow;
  int i = 0;
  for (int dz = 0; dz < _cfg->SecZ(); dz++) {
    for (int dr = 0; dr < _cfg->SecRadial(); dr++) {
      boost::shared_ptr<band> tmpBand (new band(i, dz, dr, _cfg->Din()/2.0+dr*_cfg->dDr(), _cfg->Din()/2.0+(dr+1)*_cfg->dDr(), _cfg->dDz()*dz, _cfg->dDz()*(dz+1)));
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
      boost::shared_ptr<band> tmpBand (new band(i, z, r, dRmin, dRmax, dZmin, dZmax));
      _bandAll.push_back(tmpBand);
      i++;
    }
  }
  
  //Put particles into band
  
  long long particleRemoved = 0;
  //Put particles
  for (int z = 0; z<_pRow->arraySize(); z++) {
    if (_pRow->particleReal(z)) {
      boost::shared_ptr<particle> partTemp = _pRow->getP(z);
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
        _pRow->disable(z);    //Disable and remove particle, if they are out of bands
        particleRemoved ++;
      }
    }
  }
  std::cerr<<particleRemoved<<" particles removed"<<std::endl;
 
 
 //Put forces into band
  
  long long forceRemoved = 0;
  //Put forces
  for (int z = 0; z<_fRow->arraySize(); z++) {
    boost::shared_ptr<force> forceTemp = _fRow->getF(z);
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
      _fRow->disable(z);    //Disable and remove forces, if they are out of bands
      forceRemoved ++;
    }
  }
  
  std::cerr<<particleRemoved<<" forces removed"<<std::endl;
  
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
  for(unsigned int i=0; i<_bandAll.size(); i++) {
    _bandAll[i]->calculateValues();
  }
};


void band::calculateValues () {
  _tau = 0.0;
  _p = 0.0;
  for(unsigned long long f=0; f<_allForces.size(); f++) {
    _tau += _allForces[f]->Tau();
    _p += _allForces[f]->Press();
  }
  _tauavg = _tau/_vol;
  _pavg = _p/_vol;
  
  unsigned long long i = 0;
  double angVelTmp = 0.0;
  for(unsigned long long p=0; p<_allPart.size(); p++) {
    if (not(_allPart[p]->disabled())) {
      angVelTmp += _allPart[p]->realAngular();
      i++;
    }
  }
  _vavg = angVelTmp / i;
};

