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

#include "snapshotRow.h"
#include <iostream>

std::shared_ptr<snapshot> snapshotRow::getSnapshot(unsigned int i) {
  if (i<_snapshotRow.size()) { 
    return _snapshotRow[i];
  } else {
    std::cerr<<"Illegal access to snapshotRow."<<std::endl;
    exit (EXIT_FAILURE);
  }
};

bool snp_ptr_less( std::shared_ptr<snapshot> lhs, std::shared_ptr<snapshot> rhs ) {
    return lhs->timeStep() < rhs->timeStep();
};

void snapshotRow::sortRow() {
  std::sort( _snapshotRow.begin(), _snapshotRow.end(), & snp_ptr_less );
  for(unsigned int i=0; i<_snapshotRow.size(); i++) {
    _snapshotRow[i]->id(i);
  }
};

unsigned long long snapshotRow::timeStepMin() const {
  return _snapshotRow[0]->timeStep();
};

unsigned long long snapshotRow::timeStepMax() const {
  return _snapshotRow[_snapshotRow.size()-1]->timeStep();
};

unsigned long long snapshotRow::timeStepAvg() const {
  if (_snapshotRow.size()==1) {
    return timeStepMin();
  } else {
    return (timeStepMin() + timeStepMax())/2.0;
  }
};

double snapshotRow::timeMin() const {
  return _snapshotRow[0]->time();
};

double snapshotRow::timeMax() const {
  return _snapshotRow[_snapshotRow.size()-1]->time();
};

double snapshotRow::timeAvg() const {
  if (_snapshotRow.size()==1) {
    return timeMin();
  } else {
    return (timeMin() + timeMax())/2.0;
  }
};

