#include<iostream>
#include<ctime>
#include<vector>

#include"../../src/BinaryIO.h"
#include"../../src/Analyze.h"

#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("binary,b", boost::program_options::value<string> (), "specify binary input file")
        ("restart,r", boost::program_options::value<string> (), "specify restart file")
        ("analyze,a", "analyze for Mo, Eo and Re")
        ("tecplot,t", boost::program_options::value<string> (), "specify output file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help")){
        cout << desc << endl;
        return 1;
    }

    // create a Lattice
    Lattice meins;
    Preprocess prepro;
    Timetrack timetrack;
    string filename_in = "restart.bin";

    string filename_out; // = "alternative.dat";
       
    if (vm.count("restart")) {
        filename_in = vm["restart"].as<string>();
        cout << "Restart file is: " << filename_in << ".\n" << endl ;
        bool tmpB = read_restart_file(meins, prepro, timetrack, filename_in);
        if(tmpB == false){
            cout << "failed to read input file\nno output written" << endl;
            return 1;
        }
    }

    if (vm.count("binary")) {
        filename_in = vm["binary"].as<string>();
        cout << "Binary input file is: " << filename_in << ".\n" << endl ;
        bool tmpB = read_binary(meins, filename_in);
        if(tmpB == false){
            cout << "failed to read input file\nno output written" << endl;
            return 1;
        }
    }

    if (!(vm.count("binary") || vm.count("restart"))) { 
        cout << "no input file specified\nno output written" << endl;
        return 1;
    }

    if(vm.count("analyze")){
        double resolution = prepro.getResolution();
        ParamSet params = meins.getParams();
        const double Eo = getEotvos(params, resolution);
        const double Mo = getMorton(params);
        const double Re = getReynolds(meins, resolution);
        cout << "\nEo = " << Eo << "\nMo = " << Mo << "\nRe = " << Re << endl;
    }

    if (vm.count("tecplot")) {
        filename_out = vm["tecplot"].as<string>();
        write_techplot_output_alternative(meins, filename_out);
        cout<<"\noutput written to "<< filename_out <<endl;
    }
    else {
        cout << "no output file specified\nno output written" << endl;
    }

    
    

    return 0;
}