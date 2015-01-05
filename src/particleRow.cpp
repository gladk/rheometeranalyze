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

#include "particleRow.h"
#include <iostream>

particleRow::particleRow(long long partN ) {
  _tmpP = std::make_shared<particle>();
  _realPartNum = 0;
  _allPart.assign( partN, _tmpP ); 
};

void particleRow::addP(std::shared_ptr<particle> part ) {
  if (_allPart.size()<=part->id()) {
    _allPart.resize(part->id()+1, _tmpP);
  };
  if (_allPart[part->id()] == _tmpP) {
    _allPart[part->id()] = part;
    _realPartNum ++;
  } else {
    std::cerr<<"The particles have equal IDs. "<<part->id()<<" Aborting."<<std::endl;
    exit (EXIT_FAILURE);
  }
};

long long particleRow::elementsNum() {
  return _realPartNum;
};

bool particleRow::particleReal(long long id) {
  if (arraySize()<id or _allPart[id] == _tmpP or _allPart[id]->disabled()) {
    return false;
  } else {
    return true;
  }
};

void particleRow::disable(long long id) {
  if (not(_allPart[id]->disabled())) {
    _allPart[id]->disable();
    _realPartNum--;
  }
}

std::shared_ptr<particle> particleRow::getP(long long id) {
  return _allPart[id];
}

