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

class bandShearZone {
  private:
    int _id;
    std::shared_ptr<band> _bandMinPTR;
    std::shared_ptr<band> _bandMaxPTR;
  public:
    bandShearZone(std::shared_ptr<band>, std::shared_ptr<band>);
    double RPOS () {return (_bandMaxPTR->midLinedR() - _bandMinPTR->midLinedR())/2.0 + _bandMinPTR->midLinedR();};
    double ZPOS () {return _bandMaxPTR->midLinedZ();};
    double W () {return (_bandMaxPTR->midLinedR() - _bandMinPTR->midLinedR() - _bandMinPTR->dR());};
};
