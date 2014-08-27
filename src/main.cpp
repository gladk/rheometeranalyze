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

#include "main.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace keywords = boost::log::keywords;
namespace expr = boost::log::expressions;
namespace sinks = boost::log::sinks;
namespace attrs = boost::log::attributes;

bool sortFileTimeCreate(fs::path i, fs::path j) {
  return (fs::last_write_time(i) < fs::last_write_time(j));
}


void createOutputDir(const string & outputFolder, src::severity_logger< logging::trivial::severity_level > & lg) {
  if (not fs::is_directory(outputFolder)) {
    BOOST_LOG_SEV(lg, logging::trivial::info)<<"The directory " << outputFolder<< " does not exists. Creating.";
    if (fs::create_directories(outputFolder)) {
      BOOST_LOG_SEV(lg, logging::trivial::info)<<"The directory " << outputFolder<< " created.";
    }
  }  
}

using namespace std;

int main(int ac, char* av[])
{
  string configFileName;
  string particlesFileName;
  string forcesFileName;
  string outputFolder;
  
  
  // File log
  boost::shared_ptr< logging::core > core = logging::core::get();
  boost::shared_ptr< sinks::text_file_backend > backend =
    boost::make_shared< sinks::text_file_backend >(
      keywords::file_name = "rheometer.log"
    );
    
  typedef sinks::synchronous_sink< sinks::text_file_backend > sink_t;
  boost::shared_ptr< sink_t > sink(new sink_t(backend));
  
  core->add_global_attribute("TimeStamp", attrs::local_clock());
  
  sink->set_formatter (
    expr::format("[%1%] %2%")
      % expr::attr< boost::posix_time::ptime >("TimeStamp")
      % expr::xml_decor[ expr::stream << expr::smessage ]
  );
  core->add_sink(sink);
  
  // Screen log
  
  
  boost::shared_ptr< sinks::text_ostream_backend > backendScreen = boost::make_shared< sinks::text_ostream_backend >();
#if (BOOST_VERSION > 105400)
  backendScreen->add_stream(boost::shared_ptr< std::ostream >(&std::clog, boost::empty_deleter()));
#else
  backendScreen->add_stream(boost::shared_ptr< std::ostream >(&std::clog, logging::empty_deleter()));
#endif

  backendScreen->auto_flush(true);

  // Wrap it into the frontend and register in the core.
  // The backend requires synchronization in the frontend.
  typedef sinks::synchronous_sink< sinks::text_ostream_backend > sink_tScreen;
  boost::shared_ptr< sink_tScreen > sinkScreen(new sink_tScreen(backendScreen));
  core->add_sink(sinkScreen);
  
  using namespace logging::trivial;
  src::severity_logger< severity_level > lg;

  
  std::cout<<"\n\
Rheometeranalyze\n\
Copyright (C) 2013, 2014 TU Bergakademie Freiberg\nInstitute for Mechanics and Fluid Dynamics\n\
This program comes with ABSOLUTELY NO WARRANTY.\n\
";
  
  unsigned short setVtk = 0;
  bool setUtwente = false;
  bool setContact = false;
  bool setFollowContact = false;
  bool discreteAnalyze = false;
  int setWetParticle;
  int setSnapshotsNumb;
  int setBeginSnapshot;
  int setIntOri;
  try {
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help", "produce help message")
      ("config,c", po::value<string>(), "configuration file")
      ("particle,p", po::value<string>(), "particles dump file")
      ("force,f", po::value<string>(), "forces dump file")
      ("vtk,v", po::value<unsigned short>(&setVtk)->default_value(0), "create VTK-file (0 - no export by default, 1 - all data to export, 2 - only band-data to export, 3 - force (interactions) to export)")
      ("utwente,u", "create export files for UTwente, OFF by default")
      ("contact", "perform contact analyze and creating corresponding files, OFF by default")
      ("fcontact", "perform follow-contact analyze and creating corresponding files, OFF by default")
      ("wetparticle",po::value<int>(&setWetParticle)->default_value(-1), "perform analyze of wet particle and creating corresponding files, OFF by default")
      ("intori",po::value<int>(&setIntOri)->default_value(-1), "perform analyze of interaction orientations, set number of slots, OFF by default (-1)")
      ("snapshots,s",po::value<int>(&setSnapshotsNumb)->default_value(-1), "number of snapshots to analyze, ALL by default (-1)")
      ("begin,b",po::value<int>(&setBeginSnapshot)->default_value(-1), "snapshot number from which will be done an analyze, by default (-1) last snapshots will be analyzed")
      ("output,o", po::value<string>()->default_value("output"), "output folder")
      ("discrete,d", "use discrete analyze, each snapshot will be analyzed separately")
    ;
    
    po::positional_options_description p;
    p.add("config", -1);
    po::variables_map vm;        
    po::store(po::command_line_parser(ac, av).
    options(desc).positional(p).run(), vm);
    po::notify(vm);  
    
    if (vm.count("help")) {
      cout << desc ;
      return 0;
    }
    
    if (setVtk==1) {
      BOOST_LOG_SEV(lg, info) << "VTK-file will be created, all data to export" ;
    } else if (setVtk==2) {
      BOOST_LOG_SEV(lg, info) << "VTK-file will be created, only band-data to export" ;
    } else if (setVtk==3) {
      BOOST_LOG_SEV(lg, info) << "VTK-file will be created, force (interactions) will be exported" ;
    } else if (setVtk==0) {
      BOOST_LOG_SEV(lg, info) << "VTK-file will NOT be created" ;
    } else {
      BOOST_LOG_SEV(lg, fatal) << "-v parameter can only except values 0<v<=2.\n"; 
      exit (EXIT_FAILURE);
    }
    
    if (vm.count("utwente")) {
      BOOST_LOG_SEV(lg, info) << "UTwente-files will be created" ;
      setUtwente = true;
    } else {
      BOOST_LOG_SEV(lg, info) << "UTwente-files will NOT be created" ;
    }

    if (vm.count("contact")) {
      BOOST_LOG_SEV(lg, info) << "Contact-analyze will be performed" ;
      setContact = true;
    } else {
      BOOST_LOG_SEV(lg, info) << "Contact-analyze will NOT be performed" ;
    }

    if (vm.count("fcontact")) {
      BOOST_LOG_SEV(lg, info) << "Follow-contact-analyze will be performed" ;
      setFollowContact = true;
    } else {
      BOOST_LOG_SEV(lg, info) << "Follow-contact-analyze will NOT be performed" ;
    }
    
    if (vm.count("config")) {
      BOOST_LOG_SEV(lg, info) << "config file is: " << vm["config"].as<string>() ;
    } else {
      BOOST_LOG_SEV(lg, fatal) << "config file is required, use `-c` option for that or `--help`.\n"; 
      exit (EXIT_FAILURE);
    }
    configFileName = vm["config"].as<string>();
    
    if (vm.count("output")) {
      BOOST_LOG_SEV(lg, info) << "output base folder: " << vm["output"].as<string>() ;
    }
    outputFolder = vm["output"].as<string>();

    if (vm.count("particle")) {
      BOOST_LOG_SEV(lg, info) << "particles dump-file is: " << vm["particle"].as<string>() ;
    } else {
      BOOST_LOG_SEV(lg, fatal) << "particles dump-file is required, use `-p` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    
    particlesFileName = vm["particle"].as<string>();
    
    if (vm.count("force")){
      BOOST_LOG_SEV(lg, info) << "forces dump-file is: "  << vm["force"].as<string>() ;
    } else {
      BOOST_LOG_SEV(lg, fatal) << "force dump-file is required, use `-f` option for that or `--help` for help.\n"; 
      exit (EXIT_FAILURE);
    }
    forcesFileName = vm["force"].as<string>();
    
        
    if (vm.count("discrete")) {
      BOOST_LOG_SEV(lg, info) << "Discrete analyze will be done" ;
      discreteAnalyze = true;
    }
    
  }
  catch(exception& e) {
      BOOST_LOG_SEV(lg, fatal) << "error: " << e.what() ;
      exit (EXIT_FAILURE);
  }
  catch(...) {
      BOOST_LOG_SEV(lg, fatal) << "Exception of unknown type!\n";
  }
  
  if (not(fs::is_regular_file(configFileName))) {
    fs::path p = configFileName;
    BOOST_LOG_SEV(lg, fatal)<<"The file "<<configFileName<<" does not exist. Exiting.";
    exit (EXIT_FAILURE);
  }
  
  #ifdef ALGLIB
    BOOST_LOG_SEV(lg, info)<<"ALGLIB Library is found and export of shearbands will be produced \n";
  #else
    BOOST_LOG_SEV(lg, warning)<<"ALGLIB Library is NOT found and export of shearbands will NOT be produced \n";
  #endif
  
  //=====================================================
  std::vector< fs::path > filesParticle;
  
  fs::path particle_path( particlesFileName );
  fs::path particle_dir = particle_path.parent_path();
  fs::path particle_filesmask = particle_path.filename();
  
  if (not(fs::is_directory(particle_dir))) {
    BOOST_LOG_SEV(lg, fatal)<<"The Directory "<<particle_dir.string()<<" does not exists. Exiting.";
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
    BOOST_LOG_SEV(lg, fatal)<<"The file "<<particlesFileName<<" does not exists. Exiting.";
    exit (EXIT_FAILURE);
  }
  
  sort(filesParticle.begin(), filesParticle.end(), sortFileTimeCreate); 
  
  //=====================================================
  
  std::vector< fs::path > filesForces;
    
  fs::path force_path( forcesFileName );
  fs::path force_dir = force_path.parent_path();
  fs::path force_filesmask = force_path.filename();
  
  if (not(fs::is_directory(force_dir))) {
    BOOST_LOG_SEV(lg, fatal)<<"The Directory "<<force_dir.string()<<" does not exists. Exiting.";
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
    BOOST_LOG_SEV(lg, fatal)<<"The file "<<forcesFileName<<" does not exists. Exiting.";
    exit (EXIT_FAILURE);
  } 
  
  sort(filesForces.begin(), filesForces.end(), sortFileTimeCreate); 
  
  //=====================================================
  
  if (filesParticle.size() != filesForces.size()) {
    BOOST_LOG_SEV(lg, fatal)<<"The number of force ("<<filesForces.size()<<") and particle ("<<filesParticle.size()<<") files is not the same! Exiting.";
    exit (EXIT_FAILURE);
  } else {
    BOOST_LOG_SEV(lg, info)<<"Number of particle files is "<< filesParticle.size()   ;
    BOOST_LOG_SEV(lg, info)<<"Number of force files is "<< filesForces.size()   ;
    int snapshotsNumbTemp = setSnapshotsNumb;
    int beginSnapshotTemp = setBeginSnapshot;
    
    if (setSnapshotsNumb > 0){
      if ((unsigned) setSnapshotsNumb <= filesParticle.size()) {
        if (setBeginSnapshot>=0) {
          if ((unsigned) (setBeginSnapshot+setSnapshotsNumb-1) > filesParticle.size()) {
            BOOST_LOG_SEV(lg, fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb<<", but starting from "
                      << setBeginSnapshot <<", it is only possible to have only "<< 
                      filesForces.size()-setBeginSnapshot + 1<<
                      "! Exiting.";
            exit (EXIT_FAILURE);
          }
        }
      } else if ((unsigned) setSnapshotsNumb > filesParticle.size()){
        BOOST_LOG_SEV(lg, fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb<<", but its total number is "
                      << filesForces.size() <<
                      "! Exiting.";
            exit (EXIT_FAILURE);
      }
    } else {
      if ( (setBeginSnapshot>0) and ((unsigned) setBeginSnapshot > filesParticle.size()) ) {
        BOOST_LOG_SEV(lg, fatal)<<"Requested starting snapshot is "<<setBeginSnapshot<<", but  but its total number is "
                 << filesParticle.size() <<
                 "! Exiting.";
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
      BOOST_LOG_SEV(lg, fatal)<<"Requested number of analyzed snapshots is "<<setSnapshotsNumb
               <<"! Exiting.";
      exit (EXIT_FAILURE);
    }
    
    if (beginSnapshotTemp<0) {
      beginSnapshotTemp = filesParticle.size() - snapshotsNumbTemp + 1;
    }
    
    if (filesParticle.size() > (unsigned) snapshotsNumbTemp) {
      BOOST_LOG_SEV(lg, info)<<"Reducing the number of files from "<< filesParticle.size() <<" to " << snapshotsNumbTemp ;
    }
    BOOST_LOG_SEV(lg, info)<<"Starting analyze from snapshot "<< beginSnapshotTemp ;
    filesParticle.erase(filesParticle.begin(), filesParticle.begin() + beginSnapshotTemp - 1);
    filesParticle.erase(filesParticle.begin() + snapshotsNumbTemp, filesParticle.end());
    
    filesForces.erase(filesForces.begin(), filesForces.begin() + beginSnapshotTemp - 1);
    filesForces.erase(filesForces.begin() + snapshotsNumbTemp, filesForces.end());
  }
  //=====================================================
  
  std::shared_ptr<configopt> configParams (new configopt(configFileName));
  
  if (setVtk > 0) configParams->setVtk(setVtk);
  if (setContact) configParams->setContact();
  if (setFollowContact) configParams->setFollowContact();
  if (setWetParticle) configParams->setWetParticle();
  if (setUtwente) configParams->setUtwente();
  if (setIntOri>0) configParams->setIntOri(setIntOri);
  if (setWetParticle>0) configParams->setWetParticle(setWetParticle);
  
  if (discreteAnalyze) {
    for(unsigned int i=0; i<filesParticle.size(); i++) {
      const string outputFolderNew = outputFolder + '/' + filesParticle[i].stem().string();
      createOutputDir(outputFolderNew, lg);
      std::shared_ptr<snapshotRow> snapshots (new snapshotRow());
      std::shared_ptr<snapshot> snapshotTmp (new snapshot(filesParticle[i], filesForces[i], 0));
      snapshots->addSnapshot(snapshotTmp);
      
      configParams->setSnapshot(snapshots);
      configParams->FOutput(outputFolderNew);
      
      std::shared_ptr<rheometer> curRheom (new rheometer(configParams));
    }
  } else {
    std::shared_ptr<snapshotRow> snapshots (new snapshotRow());
    createOutputDir(outputFolder, lg);
    for(unsigned int i=0; i<filesParticle.size(); i++) {
      std::shared_ptr<snapshot> snapshotTmp (new snapshot(filesParticle[i], filesForces[i], 0));
      snapshots->addSnapshot(snapshotTmp);
    }
    
    configParams->setSnapshot(snapshots);
    configParams->FOutput(outputFolder);
    
    std::shared_ptr<rheometer> curRheom (new rheometer(configParams));
  }
  fs::rename("rheometer.log", outputFolder + "/rheometer.log");
  return 0;
}
