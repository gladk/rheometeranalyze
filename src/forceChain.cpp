/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013, 2014 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, 2014 Anton Gladky <gladky.anton@gmail.com>

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

#include "forceChain.h"
#include "math_custom.h"
#include <iostream>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/foreach.hpp>

forceChain::forceChain (const std::vector <std::shared_ptr<particle> > & p, const std::vector <std::shared_ptr<force> > & f) {
  
  // Find AVG sigma3
  using namespace boost::accumulators;
  accumulator_set<double, stats<tag::mean > > acc_sigma3;
  
  BOOST_FOREACH(std::shared_ptr <particle> i,  p) {
    if (not(i->disabled())) acc_sigma3(fabs(i->stressSigma3()));
  }
  const double sigma3AVG = mean(acc_sigma3);
  
  // Mark highstressed particles
  
  BOOST_FOREACH(std::shared_ptr <particle> i,  p) {
    if (not(i->disabled()) and i->stressSigma3()>=sigma3AVG) i->highStress(1);
  }
  
  // Create pool of stressed particles, which are having >= 3 contacts with other 
  // highstressed particles
  
  
  BOOST_FOREACH(std::shared_ptr <particle> i,  p) {
    // std::cerr<<i->highStressedContacts()<<std::endl;
    if (i->highStress()>0 and i->highStressedContacts() >= 3) {
      i->highStress(i->highStressedContacts());
      _highStressPart.push_back(i);
    }
  }
  
  
}
