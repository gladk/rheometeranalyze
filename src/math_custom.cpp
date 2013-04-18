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
    return Eigen::Vector3f(rho, z, 0.0);
  } else if (x >= 0) {
    return Eigen::Vector3f(rho, z, (asin(y/rho)));
  } else if (x < 0) {
    return Eigen::Vector3f(rho, z, (-asin(y/rho) + M_PI));
  } else {
    exit (EXIT_FAILURE);
    return Eigen::Vector3f::Zero();
  }
};

Eigen::Matrix3f cart_to_cyl(Eigen::Matrix3f cart, double phi) {
  double const& axx = cart(0);
  double const& axy = cart(1);
  double const& axz = cart(2);
  double const& ayx = cart(3);
  double const& ayy = cart(4);
  double const& ayz = cart(5);
  double const& azx = cart(6);
  double const& azy = cart(7);
  double const& azz = cart(8);
  
  double const& arr = axx * cos(phi) * cos(phi) + (axy + ayz ) * cos(phi) * sin (phi) + ayy * sin(phi) * sin (phi);
  double const& arf = axy * cos(phi) * cos(phi) - (axx - ayy ) * cos(phi) * sin (phi) - ayx * sin(phi) * sin (phi);
  double const& arz = axz * cos(phi) + ayz * sin(phi);
  double const& afr = ayx * cos(phi) * cos(phi) - (axx - ayy ) * cos(phi) * sin (phi) - axy * sin(phi) * sin (phi);
  double const& aff = ayy * cos(phi) * cos(phi) - (axy + ayx ) * cos(phi) * sin (phi) + axx * sin(phi) * sin (phi);
  double const& afz = ayz * cos(phi) - axz * sin (phi);
  double const& azr = azx * cos(phi) + azy * sin (phi);
  double const& azf = azy * cos(phi) + azx * sin (phi);
  
  Eigen::Matrix3f tmpMatrix;
                  tmpMatrix << arr, arz, arf,
                               azr, azz, azf,
                               afr, afz, aff;
  return tmpMatrix;
}
