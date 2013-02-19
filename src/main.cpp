#include "main.h"

namespace po = boost::program_options;
namespace fs = boost::filesystem;
using namespace std;

int main(int ac, char* av[])
{
  string configFileName;
  string particlesFileName;
  string forcesFileName;
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

    if (vm.count("config"))
    {
        cout << "config file is: " 
              << vm["config"].as<string>() << "\n";
    }
    configFileName = vm["config"].as<string>();

    if (vm.count("particle"))
    {
        cout << "particles dump-file is: " 
              << vm["particle"].as<string>() << "\n";
    }
    
    particlesFileName = vm["particle"].as<string>();
    
    if (vm.count("force"))
    {
        cout << "forces dump-file is: " 
              << vm["force"].as<string>() << "\n";
    }
    forcesFileName = vm["force"].as<string>();

  }
  catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      return 1;
  }
  catch(...) {
      cerr << "Exception of unknown type!\n";
  }
  
  
  
  if (not(fs::is_regular_file(configFileName))) {
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
