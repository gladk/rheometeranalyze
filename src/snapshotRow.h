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
#include "snapshot.h"

class snapshotRow {
  private:
    std::vector <std::shared_ptr<snapshot>> _snapshotRow;
  public:
    void addSnapshot(std::shared_ptr<snapshot> snapshotTmp){_snapshotRow.push_back(snapshotTmp);};
    unsigned int size() const {  return _snapshotRow.size();};
    std::shared_ptr<snapshot> getSnapshot(unsigned int i);
    void sortRow();
    unsigned long long timeStepMin() const;
    unsigned long long timeStepMax() const;
    unsigned long long timeStepAvg() const;
    double timeMin() const;
    double timeMax() const;
    double timeAvg() const;
    void showTimes() const;
    void clear();
};
