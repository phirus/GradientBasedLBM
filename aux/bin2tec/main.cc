#include<iostream>
#include<ctime>
#include<vector>

#include"../../src/BinaryIO.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("input,i", boost::program_options::value<string> (), "specify restart file")
        ("output,o", boost::program_options::value<string> (), "specify output file")
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

    string filename_out = "alternative.dat";
       
    if (vm.count("input")) {
        filename_in = vm["input"].as<string>();
        cout << "Restart file is: " << vm["input"].as<string>() << ".\n" << endl ;
    }

    if (vm.count("output")) {
        filename_out = vm["output"].as<string>();
    }

    bool tmpB = read_restart_file(meins, prepro, timetrack, filename_in);
    if(tmpB == false){
        cout << "failed to read input file\nno output written" << endl;
        return 1;
    }

    write_techplot_output_alternative(meins, filename_out);
    cout<<"\noutput written to "<< filename_out <<endl;

    return 0;
}