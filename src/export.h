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

#pragma once

#include <Eigen/Dense>
#include <memory>

#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkZLibDataCompressor.h>
#include <vtkLine.h>

#include "band.h"
#include "bandRow.h"
#include "forceRow.h"
#include "forceChain.h"
#include "math_custom.h"

struct contactFollow {
  std::shared_ptr<force> _f;     // Pointer into force
  std::shared_ptr<snapshot> _sn; // Pointer into snapshot
    contactFollow(std::shared_ptr<force> f, std::shared_ptr<snapshot> sn) {_f = f; _sn = sn;};
  unsigned long long P1_id() {return _f->pid1();};
  unsigned long long P2_id() {return _f->pid2();};
  unsigned long long timeStep() {return _sn->timeStep();};
  double deltaN() {return _f->deltaN();};
  double volWater() {return _f->volWater();};
  double distCurr() {return _f->distCurr();};
  double distCrit() {return _f->distCrit();};
  Eigen::Vector3d f() {return _f->val();};
};

bool sortContactFollow(std::shared_ptr<contactFollow> i, std::shared_ptr<contactFollow> j);

class exportclass {
  private:
    std::shared_ptr <configopt> _cfg;
    std::shared_ptr <bandRow> _bandRow;
    std::vector<std::shared_ptr <forceRow>> _forceAll;
    
  public:
    exportclass(const std::shared_ptr<configopt>, const std::shared_ptr <bandRow>, const std::vector<std::shared_ptr <forceRow>>);
    void VTK();
    void gnuplotSchearRate();
    void gnuplotContactAnalyze(int bins);
    void gnuplotContactNumberAnalyze();
    void gnuplotContactWet();
    void gnuplotContactFollow();
    void gnuplotWetParticles(int bins);
    void Utwente();
    void intOri();
    void torque();
    void forceChain();
};
