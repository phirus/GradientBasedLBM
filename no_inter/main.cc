#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../src/2D/Lattice2D_no_inter.h"
#include"../src/2D/BinaryIO2D.h"
#include"../src/2D/Analyze2D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);
void initialSetUp_Diffusion(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);
void initialSetUp_Membrane(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);
void initialSetUp_Chamber(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);
void initializeShearfFlow(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("boundary,b", boost::program_options::value<string> (), "specify boundary input file")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ("override,o", boost::program_options::value<string> (), "specify parameter file to bypass preprocess routine")
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
    ParamSet params = prepro.getParamSet();
    Boundaries boundaries = read_boundaries_file("BoundaryInput");

    Preprocess prepro_input;
    Timetrack timetrack_input;    
    ParamSet params_input;

    if (vm.count("preprocess")) 
    {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro_input = read_preprocess_file(vm["preprocess"].as<string>());
        timetrack_input = read_timetrack_file(vm["preprocess"].as<string>());
        params_input = prepro_input.getParamSet();
        
        prepro = prepro_input;
        timetrack = timetrack_input;
        params = params_input;
    }

        if (vm.count("boundary")) 
    {
        cout << "boundary file is: " << vm["boundary"].as<string>() << ".\n" << endl ;
        boundaries = read_boundaries_file(vm["boundary"].as<string>());
    }

    if (vm.count("override")) 
    {
        cout << "override preprocess file with: " << vm["override"].as<string>() << ".\n" << endl ;
        params_input = read_paramset_file(vm["override"].as<string>());
        params = params_input;
    }

    // create a Lattice   
    int ymax = prepro.getYCells();
    int xmax = prepro.getXCells();
    Lattice2D_no_inter meins(xmax,ymax);

//    initialSetUp(meins, prepro, boundaries, xmax, ymax, params);
    initialSetUp_Membrane(meins, prepro, boundaries, xmax, ymax, params);
    write_vtk_output2D(meins, 0);
    

    const int outputInterval = timetrack.getOutputInt();


    int i;

    while (timetrack.proceed() == true)
    {
        meins.collideAll(numOfCPUs);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%10 == 0) 
        {
            cout << i <<endl;
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        }        
    }


    return 0;
}


void initialSetUp(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    const Cell2D air(0,1.0,false);
    const Cell2D liquid(1.0,0,false);

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax * 0.5;
    const int ym1 = ymax * 0.5;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.setBoundaries(bound);
    meins.buildWalls();
    meins.overallRho();
  
    //meins.equilibriumIni();
    cout<<"Initialisierung beendet\n\n"<<endl;
}

void initialSetUp_Diffusion(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    const Cell2D air(0,1.0,false);
    const Cell2D liquid(1.0,0,false);
    const Cell2D wall(0,0,true);

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax * 0.5;
    const int ym1 = ymax / 4;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);

            if (j == ymax/2 && i%3 != 0) meins.setCell(i,j,wall);

        }
    }

    meins.setBoundaries(bound);
    meins.buildWalls();
    meins.overallRho();
  
    //meins.equilibriumIni();
    cout<<"Initialisierung beendet\n\n"<<endl;
}

void initialSetUp_Membrane(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    const Cell2D air(0,1.01,false);
    const Cell2D liquid(1.0,0,false);
    const Cell2D wall(0,0,true);

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax / 2;
    const int ym1 = ymax / 2;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            // if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);

            // if (j == ymax/2 && i%4 != 0) meins.setCell(i,j,wall);
            // if (j == ymax/2 && i!= xmax/2) meins.setCell(i,j,wall);

        }
    }

    meins.setBoundaries(bound);
    meins.buildWalls();
    meins.overallRho();
  
    //meins.equilibriumIni();
    cout<<"Initialisierung beendet\n\n"<<endl;
}

void initialSetUp_Chamber(Lattice2D_no_inter& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    const Cell2D air(0,1.0,false);
    const Cell2D liquid(1.0,0,false);
    const Cell2D wall(0,0,true);

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax / 2;
    const int ym1 = ymax / 2;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++)
        {
            meins.setCell(i,j,liquid);
            if(i == xmax / 5 && j < ymax/5 * 4) meins.setCell(i,j,wall);
            if(i == xmax / 5 * 4 && j < ymax/5 * 4) meins.setCell(i,j,wall);
            if(j == 0 && i < xmax / 5) meins.setCell(i,j,wall);
            if(j == 0 && i > xmax / 5 * 4) meins.setCell(i,j,wall);
        }
    }

    meins.setBoundaries(bound);
    meins.buildWalls();
    meins.overallRho();
  
    //meins.equilibriumIni();
    cout<<"Initialisierung beendet\n\n"<<endl;
}