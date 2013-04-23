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

#include <Eigen/Dense>

// http://en.wikipedia.org/wiki/Cylindrical_coordinate_system#Cartesian_coordinates
Eigen::Vector3d cart_to_cyl(Eigen::Vector3d cart);                //Returns (rho, z, phi)

/*
  @book{schade2007strömungslehre,
  title={Str{\"o}mungslehre},
  author={Schade, H. and Kunz, E.},
  isbn={9783110189728},
  series={De Gruyter Lehrbuch},
  url={http://books.google.de/books?id=CbboC1qhMsUC},
  year={2007},
  publisher={Walter De Gruyter Incorporated}
  }
* Seite 517
*/ 
Eigen::Matrix3d cart_to_cyl(Eigen::Matrix3d& cart, double& phi);    //Returns (rho-rho  rho-z  rho-phi)
                                                                  //        (z-rho    z-z    z-phi  )
                                                                  //        (phi-rho  phi-z  phi-phi)

// Returns cartesians coordinates of every unit-vector of cylindrical coordinates. 
// http://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates
Eigen::Matrix3d get_axes_coord(Eigen::Vector3d& X, Eigen::Quaternion<double>& rotateCC);  //Returns [drX, drY, drZ]
                                                                                               // [dzX, dzY, dzZ]
                                                                                               // [dfX, dfY, dfZ]
// Returns vectors, which are projections on cylindrical coordinates.
Eigen::Vector3d get_cyl_rotated_vector(Eigen::Vector3d& X, Eigen::Matrix3d& rotateMatrix);  //Returns [dr, dz, drf]
