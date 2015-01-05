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

#include "forceRow.h"
#include "particleRow.h"

class bandRowBase {
  protected:
    std::vector<std::shared_ptr<band> > _bandAll;
    std::vector<Eigen::Vector3d> _shearBands;       //[Rz, W, H]
    std::shared_ptr<configopt> _cfg;
    std::vector<double> _omega0;
    double _omega0AVG=0.0;
  public:
    int getBandR(double);
    int getBandZ(double);
    int getBandF(double);
    unsigned int shearBandSize() {return _shearBands.size();};
    Eigen::Vector3d shearBand(unsigned int i) {return _shearBands[i];};
    std::shared_ptr<band> getBand(unsigned int id) {return _bandAll[id];}
    std::shared_ptr<band> getBand(unsigned int idR, unsigned int idZ) {return _bandAll[idZ*_cfg->SecRadial() + idR];}
    unsigned int size() {return _bandAll.size();}
    double totalVolume();
    double shearBandVolume();
    double omega0AVG() const;
};

class bandRow : public bandRowBase{
  private:
    std::vector<std::shared_ptr<particleRow>> _pRow;
    std::vector<std::shared_ptr<forceRow>> _fRow;
  public:
    bandRow     (std::shared_ptr<configopt> cfg, std::vector<std::shared_ptr<particleRow>> pRow, std::vector<std::shared_ptr<forceRow>> fRow);
    void fillBands();
    void calculateValues();
};
