/*
    This file is part of Rheometeranalyze.
    Rheometeranalyze is the programm to analyze the DEM-simulations of rheometer.
    Copyright (C) 2013 TU Bergakademie Freiberg, Institute for Mechanics and Fluid Dynamics

    Author: 2013, Anton Gladky <gladky.anton@gmail.com>

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

#include "main.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

using namespace std;

int main(int ac, char* av[])
{
  string configFileName;
  string particlesFileName;
  string forcesFileName;
  string outputFolder;
  
  std::cout<<"\n\
Rheometeranalyze\n\
Copyright (C) 2013 TU Bergakademie Freiberg\nInstitute for Mechanics and Fluid Dynamics\n\
This program comes with ABSOLUTELY NO WARRANTY.\n\
"<<std::endl;
  
  bool setVtk = false;
  bool setUtwente = false;
  int setSnapshotsNumb;
  int setBeginSnapshot;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("config,c", po::value<string>(), "configuration file")
      ("particle,p", po::value<string>(), "particles dump file")
      ("force,f", po::value<string>(), "forces dump file")
      ("vtk,v", "create VTK-file, OFF by default")
      ("utwente,u", "create export files for UTwente, OFF by default")
      ("snapshots,s",po::value<int>(&setSnapshotsNumb)->default_value(-1), "number of snapshots to analyze, ALL by default (-1)")
      ("begin,b",po::value<int>(&setBeginSnapshot)->default_value(-1), "snapshot number from which will be done an analyze, by default (-1) last snapshots will be analyzed")
      ("output,o", po::value<string>()->default_value("output"), "output folder")
    ;
    
    po::positional_options_description p;
    p.add("config", -1);
    po::variables_map vm;        
    po::store(po::command_line_parser(ac, av).
    options(desc).positional(p).run(), vm);
    po::notify(vm);  
    
    if (vm.count("help")) {
      cout << desc << std::endl;
      return 0;
    }
    
    if (vm.count("vtk")) {
      BOOST_LOG_TRIVIAL(trace) << "VTK-file will be created" << std::endl;
      setVtk = true;
    } else {
      BOOST_LOG_TRIVIAL(trace) << "VTK-file will NOT be created" << std::endl;
    }
    
    if (vm.count("utwente")) {
      BOOST_LOG_TRIVIAL(trace) << "UTwente-files will be created" << std::endl;
      setUtwente = true;
    } else {
      BOOST_LOG_TRIVIAL(trace) << "UTwente-files will NOT be created" << std::endl;
    }

    if (vm.count("config")) {
      BOOST_LOG_TRIVIAL(trace) << "config file is: " << vm["config"].as<string>() << std::endl;
    } else {
      BOOST_LOG_TRIVIAL(fatal) << "config file is required, use `-c` option for that or `--help`.\n"; 
      exit (EXIT_FAILURE);
    }
    configFileName = vm["config"].as<string>();
    
    if (vm.count("output")) {
      BOOST_LOG_TRIVIAL(trace) << "output folder: " << vm["output"].as<string>() << std::endl;
    }
    outputFolder = vm["output"].as<string>();

    if (vm.count("particle")) {
      BOOST_LOG_TRIVIAL(trace) << "particles dump-file is: " << vm["particle"].as<string>() << std::endl;
    } else {
      BOOST_LOG_TRIVIAL(fatal) << "particles dump-file is required, use `-p` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    
    particlesFileName = vm["particle"].as<string>();
    
    if (vm.count("force")){
      BOOST_LOG_TRIVIAL(trace) << "forces dump-file is: "  << vm["force"].as<string>() << std::endl;
    } else {
      BOOST_LOG_TRIVIAL(fatal) << "force dump-file is required, use `-f` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    forcesFileName = vm["force"].as<string>();

  }
  catch(exception& e) {
      BOOST_LOG_TRIVIAL(fatal) << "error: " << e.what() << std::endl;
      exit (EXIT_FAILURE);
  }
  catch(...) {
      BOOST_LOG_TRIVIAL(fatal) << "Exception of unknown type!\n";
  }
  
  if (not(fs::is_regular_file(configFileName))) {
    fs::path p = configFileName;
    BOOST_LOG_TRIVIAL(fatal)<<"The file "<<configFileName<<" does not exist. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  #ifdef ALGLIB
    BOOST_LOG_TRIVIAL(trace)<<"\nALGLIB Library is found and export of shearbands will be produced \n"<<std::endl;
  #else
    BOOST_LOG_TRIVIAL(warning)<<"\nALGLIB Library is NOT found and export of shearbands will NOT be produced \n"<<std::endl;
  #endif
  
  
  //=====================================================
  std::vector< fs::path > filesParticle;
  
  fs::path particle_path( particlesFileName );
  fs::path particle_dir = particle_path.parent_path();
  fs::path particle_filesmask = particle_path.filename();
  
  if (not(fs::is_directory(particle_dir))) {
    BOOST_LOG_TRIVIAL(fatal)<<"The Directory "<<particle_dir.string()<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  } else {
    fs::directory_iterator end_itr; // Default ctor yields past-the-end
    const boost::regex my_filter( particle_filesmask.string() );
    for( boost::filesystem::directory_iterator i( particle_dir.string() ); i != end_itr; ++i )
    {
      // Skip if not a file
      if( !boost::filesystem::is_regular_file( i->status() ) ) continue;
      boost::smatch what;
      // Skip if no match
      if( !boost::regex_match(  i->path().filename().string(), what, my_filter ) ) continue;
      // File matches, store it
      if (fs::is_regular_file(i->path())) {
        filesParticle.push_back(i->path());
      }
    }
  }
  
  if (filesParticle.size()<1) {
    BOOST_LOG_TRIVIAL(fatal)<<"The file "<<particlesFileName<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  }
  
  sort(filesParticle.begin(), filesParticle.end()); 
  
  //=====================================================
  
  std::vector< fs::path > filesForces;
    
  fs::path force_path( forcesFileName );
  fs::path force_dir = force_path.parent_path();
  fs::path force_filesmask = force_path.filename();
  
  if (not(fs::is_directory(force_dir))) {
    BOOST_LOG_TRIVIAL(fatal)<<"The Directory "<<force_dir.string()<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  } else {
    fs::directory_iterator end_itr; // Default ctor yields past-the-end
    const boost::regex my_filter( force_filesmask.string() );
    for( boost::filesystem::directory_iterator i( force_dir.string() ); i != end_itr; ++i )
    {
      // Skip if not a file
      if( !boost::filesystem::is_regular_file( i->status() ) ) continue;
      boost::smatch what;
      // Skip if no match
      if( !boost::regex_match(  i->path().filename().string(), what, my_filter ) ) continue;
      // File matches, store it
      if (fs::is_regular_file(i->path())) {
        filesForces.push_back(i->path());
      }
    }
  }
  
  if (filesForces.size()<1) {
    BOOST_LOG_TRIVIAL(fatal)<<"The file "<<forcesFileName<<" does not exists. Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  } 
  
  sort(filesForces.begin(), filesForces.end()); 
  
  //=====================================================
  if (not fs::is_directory(outputFolder)) {
    BOOST_LOG_TRIVIAL(trace)<<"The directory " << outputFolder<< " does not exists. Creating."<<std::endl;
    if (fs::create_directory(outputFolder)) {
      BOOST_LOG_TRIVIAL(trace)<<"The directory " << outputFolder<< " created."<<std::endl;
    }
  }
  
  //=====================================================
  
  if (filesParticle.size() != filesForces.size()) {
    BOOST_LOG_TRIVIAL(fatal)<<"The number of force ("<<filesForces.size()<<") and particle ("<<filesParticle.size()<<") files is not the same! Exiting."<<std::endl;
    exit (EXIT_FAILURE);
  } else {
    BOOST_LOG_TRIVIAL(trace)<<"Number of particle files is "<< filesParticle.size()   <<std::endl;
    BOOST_LOG_TRIVIAL(trace)<<"Number of force files is "<< filesForces.size()   <<std::endl;
    int snapshotsNumbTemp = setSnapshotsNumb;
    int beginSnapshotTemp = setBeginSnapshot;
    
    if (setSnapshotsNumb > 0){
      if (setSnapshotsNumb <= filesParticle.size()) {
        if (setBeginSnapshot>=0) {
          if ((setBeginSnapshot+setSnapshotsNumb-1) > filesParticle.size()) {
            BOOST_LOG_TRIVIAL(fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb<<", but starting from "
                      << setBeginSnapshot <<", it is only possible to have only "<< 
                      filesForces.size()-setBeginSnapshot + 1<<
                      "! Exiting."<<std::endl;
            exit (EXIT_FAILURE);
          }
        }
      } else if (setSnapshotsNumb > filesParticle.size()){
        BOOST_LOG_TRIVIAL(fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb<<", but its total number is "
                      << filesForces.size() <<
                      "! Exiting."<<std::endl;
            exit (EXIT_FAILURE);
      }
    } else {
      if ( (setBeginSnapshot>0) and (setBeginSnapshot > filesParticle.size()) ) {
        BOOST_LOG_TRIVIAL(fatal)<<"Requested starting snapshot is "<<setBeginSnapshot<<", but  but its total number is "
                 << filesParticle.size() <<
                 "! Exiting."<<std::endl;
            exit (EXIT_FAILURE);
        exit (EXIT_FAILURE);
      } else if (setBeginSnapshot>0) {
        snapshotsNumbTemp = filesParticle.size() - setBeginSnapshot + 1;
      } else {
        snapshotsNumbTemp = filesParticle.size();
        beginSnapshotTemp = 1;
      }
    }
   
    if (snapshotsNumbTemp<0) {
      snapshotsNumbTemp = filesParticle.size();
    } else if (snapshotsNumbTemp==0) {
      BOOST_LOG_TRIVIAL(fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb
               <<"! Exiting."<<std::endl;
      exit (EXIT_FAILURE);
    }
    
    if (beginSnapshotTemp<0) {
      beginSnapshotTemp = filesParticle.size() - snapshotsNumbTemp + 1;
    }
    
    if (filesParticle.size() > snapshotsNumbTemp) {
      BOOST_LOG_TRIVIAL(trace)<<"Reducing the number of files from "<< filesParticle.size() <<" to " << snapshotsNumbTemp <<std::endl;
    }
    BOOST_LOG_TRIVIAL(trace)<<"Starting analyze from snapshot "<< beginSnapshotTemp <<std::endl;
    filesParticle.erase(filesParticle.begin(), filesParticle.begin() + beginSnapshotTemp - 1);
    filesParticle.erase(filesParticle.begin() + snapshotsNumbTemp, filesParticle.end());
    
    filesForces.erase(filesForces.begin(), filesForces.begin() + beginSnapshotTemp - 1);
    filesForces.erase(filesForces.begin() + snapshotsNumbTemp, filesForces.end());

    
  }
  
  //=====================================================
  
  std::shared_ptr<snapshotRow> snapshots (new snapshotRow());
  
  for(unsigned int i=0; i<filesParticle.size(); i++) {
    std::shared_ptr<snapshot> snapshotTmp (new snapshot(filesParticle[i], filesForces[i], 0));
    snapshots->addSnapshot(snapshotTmp);
  }
  
  
  std::shared_ptr<configopt> configParams (new configopt(configFileName));
  configParams->setSnapshot(snapshots);
  configParams->FOutput(outputFolder);
  
  if (setVtk) configParams->setVtk();
  if (setUtwente) configParams->setUtwente();
  
  std::shared_ptr<rheometer> curRheom (new rheometer(configParams));
  
  return 0;
}
