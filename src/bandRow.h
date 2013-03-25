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

#pragma once

#include "bandShearZone.h"
#include "forceRow.h"

class bandRow {
  private:
    std::vector<std::shared_ptr<band> > _bandAll;
    std::vector<std::shared_ptr<bandShearZone> > _bandShearZones;
    std::shared_ptr<configopt> _cfg; 
    std::vector<std::shared_ptr<particleRow>> _pRow;
    std::vector<std::shared_ptr<forceRow>> _fRow;
  public:
    bandRow (std::shared_ptr<configopt>, std::vector<std::shared_ptr <particleRow>>, std::vector<std::shared_ptr <forceRow>>);
    void fillBands();
    int getBandR(double);
    int getBandZ(double);
    unsigned int getBandShearZonesSize() {return _bandShearZones.size();};
    std::shared_ptr<bandShearZone> getBandShearZones(int id) {return _bandShearZones[id];};
    void calculateValues();
    std::shared_ptr<band> getBand(unsigned int id) {return _bandAll[id];}
    std::shared_ptr<band> getBand(unsigned int idR, unsigned int idZ) {return _bandAll[idZ*_cfg->SecRadial() + idR];}
    unsigned int size() {return _bandAll.size();}
};

