#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../src/2D/Lattice2D.h"
#include"../src/2D/BinaryIO2D.h"
#include"../src/2D/Analyze2D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("boundary,b", boost::program_options::value<string> (), "specify boundary input file")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help"))
    {
        cout << desc << endl;
        return 1;
    }

    int numOfCPUs = 1;
    if(vm.count("cpu"))
    {
        numOfCPUs = vm["cpu"].as<int>();
        cout << "number of CPUs set to "<< numOfCPUs << endl;
    }

    Preprocess prepro = read_preprocess_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file("preprocessFile");
    Boundaries boundaries = read_boundaries_file("BoundaryInput");

    if (vm.count("preprocess")) 
    {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro = read_preprocess_file(vm["preprocess"].as<string>());
        timetrack = read_timetrack_file(vm["preprocess"].as<string>());
    }

        if (vm.count("boundary")) 
    {
        cout << "boundary file is: " << vm["boundary"].as<string>() << ".\n" << endl ;
        boundaries = read_boundaries_file(vm["boundary"].as<string>());
    }

    ParamSet params = prepro.getParamSet();

    // create a Lattice   
    int ymax = 600; 
    int xmax = 600; 
    Lattice2D meins(xmax,ymax);

    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell2D air(0,rho_gas,false);
    const Cell2D liquid(rho_liquid,0,false);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = 15;
    const int xm1 = xmax / 2;
    const int ym1 = ymax / 2;

    BubbleBox2D bubblebox;
    bubblebox.setBubble(300,300);
    bubblebox.setW(100);
    bubblebox.setH(100);
    meins.setBubbleBox(bubblebox);

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }
    meins.setBoundaries(boundaries);

    for (int i = 1; i< 1001; i++)
    {
        meins.collideAll(numOfCPUs,false,false);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);
    }
    write_vtk_output2D(meins, 1);

    return 0;
}