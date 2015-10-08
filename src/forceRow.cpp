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

#include "main.h"
#include "forceRow.h"

forceRow::forceRow() {
  std::vector <std::shared_ptr<force> > _allForce;
};

void forceRow::addF(std::shared_ptr<force> forc ) {
  _allForce.push_back(forc);
};

unsigned long long forceRow::elementsNum() {
  return _allForce.size();
};

std::shared_ptr<force> forceRow::getF(unsigned long long id) {
  return _allForce[id];
};

void forceRow::clear() {
  _allForce.clear();
};
