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

#include "main.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

int main(int ac, char* av[])
{
  string configFileName;
  string particlesFileName;
  string forcesFileName;
  
  std::cerr<<"\n\
Rheometeranalyze\n\
Copyright (C) 2013 TU Bergakademie Freiberg\nInstitute for Mechanics and Fluid Dynamics\n\
This program comes with ABSOLUTELY NO WARRANTY.\n\
"<<std::endl;
  
  try {

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("config,c", po::value<string>(), "configuration file")
        ("particle,p", po::value<string>(), "particles dump file")
        ("force,f", po::value<string>(), "forces dump file")
    ;
    
    po::positional_options_description p;
    p.add("config", -1);
    po::variables_map vm;        
    po::store(po::command_line_parser(ac, av).
    options(desc).positional(p).run(), vm);
    po::notify(vm);  

    if (vm.count("help")) {
        cout << desc << "\n";
        return 0;
    }

    if (vm.count("config")) {
      cout << "config file is: " << vm["config"].as<string>() << "\n";
    } else {
      cout << "config file is required, use `-c` option for that or `--help`.\n"; 
      exit (EXIT_FAILURE);
    }
    configFileName = vm["config"].as<string>();

    if (vm.count("particle")) {
      cout << "particles dump-file is: " << vm["particle"].as<string>() << "\n";
    } else {
      cout << "particles dump-file is required, use `-p` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    
    particlesFileName = vm["particle"].as<string>();
    
    if (vm.count("force")){
      cout << "forces dump-file is: "  << vm["force"].as<string>() << "\n";
    } else {
      cout << "force dump-file is required, use `-p` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    forcesFileName = vm["force"].as<string>();

  }
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      exit (EXIT_FAILURE);
  }
  catch(...) {
      cerr << "Exception of unknown type!\n";
  }
  
  if (not(fs::is_regular_file(configFileName))) {
    fs::path p = configFileName;
    std::cerr<<"The file "<<configFileName<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  if (not(fs::is_regular_file(particlesFileName))) {
    std::cerr<<"The file "<<particlesFileName<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  if (not(fs::is_regular_file(forcesFileName))) {
    std::cerr<<"The file "<<forcesFileName<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  std::shared_ptr<configopt> configParams (new configopt(configFileName));
  std::shared_ptr<rheometer> curRheom (new rheometer(configParams, particlesFileName, forcesFileName));
  
  return 0;
}
