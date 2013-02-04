#include <boost/program_options.hpp>

#include <iostream>
#include <iterator>
#include <Eigen/Dense>

namespace po = boost::program_options;
namespace eig= Eigen;
using namespace std;

int main(int ac, char* av[])
{
    try {

        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")
            ("x-center", po::value<double>()->default_value(0.0), "set x-center coordinate")
            ("y-center", po::value<double>()->default_value(0.0), "set y-center coordinate")
            ("z-center", po::value<double>()->default_value(0.0), "set z-center coordinate")
            ("Din", po::value<double>()->default_value(0.0), "set internal diameter of the rheometer")
            ("Dout", po::value<double>()->default_value(0.0), "set external diameter of the rheometer")
            ("RadialSec", po::value<int>()->default_value(0), "set numer of sections, in which the rheometer will be divided in radial direction")
            ("ZSec", po::value<int>()->default_value(0), "set numer of sections, in which the rheometer will be divided in z-direction")
        ;

        po::variables_map vm;        
        po::store(po::parse_command_line(ac, av, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 0;
        }

        if (vm.count("x-center")) {
            cout << "x-center was set to " 
                 << vm["x-center"].as<double>() << ".\n";
        }
        
        if (vm.count("y-center")) {
            cout << "y-center was set to " 
                 << vm["y-center"].as<double>() << ".\n";
        }
        
        if (vm.count("z-center")) {
            cout << "z-center was set to " 
                 << vm["z-center"].as<double>() << ".\n";
        }

        if (vm.count("Din")) {
            cout << "Din was set to " 
                 << vm["Din"].as<double>() << ".\n";
        }
        
        if (vm.count("Dout")) {
            cout << "Dout was set to " 
                 << vm["Dout"].as<double>() << ".\n";
        }
        
        
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        return 1;
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
    }

    return 0;
}
