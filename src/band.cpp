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
};


bandRow::bandRow (boost::shared_ptr<configopt> cfg, boost::shared_ptr<particleRow> pRow){
  _cfg =  cfg;
  _pRow = pRow;
  int i = 0;
  for (int dz = 0; dz < _cfg->SecZ(); dz++) {
    for (int dr = 0; dr < _cfg->SecRadial(); dr++) {
      boost::shared_ptr<band> tmpBand (new band(i, dz, dr, _cfg->Din()/2.0+dr*_cfg->dDr(), _cfg->Din()/2.0+(dr+1)*_cfg->dDr(), _cfg->dDz()*dz, _cfg->dDz()*(dz+1)));
      _bandAll.push_back(tmpBand);
      i++;
    }
  }
  fillBands();
};
  
void bandRow::fillBands (){
  //Fill bands with particles
  Eigen::Vector3d O = _cfg->get_c();
  Eigen::Vector3d Z = _cfg->get_o();
  for (int z = 0; z<_pRow->arraySize(); z++) {
    if (_pRow->particleReal(z)) {
      boost::shared_ptr<particle> partTemp = _pRow->getP(z);
      Eigen::Vector3d OP = partTemp->c() - O;     //Vector from center to point
      Eigen::Vector3d OPV = OP.cross(Z);          //Vector, temporal
      OPV.normalize();
      Eigen::Vector3d OPV1 = Z.cross(OPV);        //Vector for projection
      double dist = OP.dot(OPV1);
      partTemp->set_dist(dist);

    }
  }
}
