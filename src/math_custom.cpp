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

#include "math_custom.h"

Eigen::Vector3d cart_to_cyl(Eigen::Vector3d cart) {
  double const& x = cart(0);
  double const& y = cart(1);
  double const& z = cart(2);
  
  double rho = sqrt(x*x + y*y);
  
  return Eigen::Vector3d(rho, z, (atan2(y,x)));

};

Eigen::Matrix3d cart_to_cyl(Eigen::Matrix3d& cart, double& phi) {
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
  
  Eigen::Matrix3d tmpMatrix;
                  tmpMatrix << arr, arz, arf,
                               azr, azz, azf,
                               afr, afz, aff;
  return tmpMatrix;
}

Eigen::Matrix3d get_axes_coord(Eigen::Vector3d& X, Eigen::Quaternion<double>& rotateCC) {
  double const& z = X(1);
  double const& phi = X(2);
  
  Eigen::Vector3d dr = Eigen::Vector3d(cos(phi), sin(phi), 0.0); dr = rotateCC*dr;
  Eigen::Vector3d df = Eigen::Vector3d(-sin(phi), cos(phi), 0.0); df = rotateCC*df;
  Eigen::Vector3d dz = Eigen::Vector3d(0.0, 0.0, z); dz = rotateCC*dz;
  dr.normalize(); dz.normalize(); df.normalize();
  
  Eigen::Matrix3d tmpMatrix;
  tmpMatrix << dr, dz, df;
  tmpMatrix.transposeInPlace();
  
  return tmpMatrix;
}

Eigen::Vector3d get_cyl_rotated_vector(Eigen::Vector3d& X, Eigen::Matrix3d& rotateMatrix) {
    
  Eigen::Matrix3d tempMatrix; tempMatrix << X, X, X;
  tempMatrix.transposeInPlace();
  tempMatrix = rotateMatrix.cwiseProduct(tempMatrix);
  
  Eigen::Vector3d tmpVector = Eigen::Vector3d(tempMatrix.row(0).sum(),
                                              tempMatrix.row(1).sum(),
                                              tempMatrix.row(2).sum());
  return tmpVector;
}

Eigen::Vector3d cart_to_sph(Eigen::Vector3d cart) {                //Input = [rho, z, phi], return [Theta, Psi, rho]
  double const& x = cart(2);
  double const& y = cart(0);
  double const& z = cart(1);
  
  double const rho = sqrt(x*x + y*y + z*z);
  double Theta = atan2(y,x);
  double const Psi = acos(z/rho);
  
  if (Theta < 0.0) {
    Theta+=2*M_PI;
  }
  return Eigen::Vector3d(Theta, Psi, rho);
}

Eigen::Vector3d sph_to_cart(Eigen::Vector3d sph) {                //Input = [Theta(0<=Theta<=2Pi), Psi(0<=Psi<=Pi), R], return = [X, Y, Z]
  double const& Theta = sph(0);
  double const& Psi   = sph(1);
  double const& rho   = sph(2);
  
  double const x = rho * sin (Psi) * cos(Theta);
  double const y = rho * sin (Psi) * sin(Theta);
  double const z = rho * cos (Psi);

  return Eigen::Vector3d(x, y, z);
}
