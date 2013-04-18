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

#include "math_custom.h"

// http://en.wikipedia.org/wiki/Cylindrical_coordinate_system#Cartesian_coordinates
Eigen::Vector3f cart_to_cyl(Eigen::Vector3f cart) {
  double const& x = cart(0);
  double const& y = cart(1);
  double const& z = cart(2);
  
  double rho = sqrt(x*x + y*y);
  
  if (x == 0.0 and y == 0.0) {
    return Eigen::Vector3f(rho, 0.0, z);
  } else if (x >= 0) {
    return Eigen::Vector3f(rho, (asin(y/rho)), z);
  } else if (x < 0) {
    return Eigen::Vector3f(rho, (-asin(y/rho) + M_PI), z);
  } else {
    exit (EXIT_FAILURE);
    return Eigen::Vector3f::Zero();
  }
};
