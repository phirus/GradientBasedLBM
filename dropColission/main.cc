#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../src/2D/Lattice2D.h"
#include"../src/2D/BinaryIO2D.h"
#include"../src/2D/Analyze2D.h"
//#include"../src/Preprocess_Drop.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice2D& meins, Preprocess_Drop prepro, Boundaries bound, int xmax, int ymax, ParamSet params, int numOfCPUs);
const Boundaries getBoundaries(double gamma);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
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

    Preprocess_Drop prepro = read_preprocess_drop_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file("preprocessFile");
    ParamSet params = prepro.getParamSet();
    
    Preprocess_Drop prepro_input;
    Timetrack timetrack_input;    
    ParamSet params_input;

    if (vm.count("preprocess")) 
    {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro_input = read_preprocess_drop_file(vm["preprocess"].as<string>());
        timetrack_input = read_timetrack_file(vm["preprocess"].as<string>());
        params_input = prepro_input.getParamSet();
        
        prepro = prepro_input;
        timetrack = timetrack_input;
        params = params_input;
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
    Lattice2D meins(xmax,ymax);

    Boundaries bound ;//= getBoundaries(prepro.getGamma());

    //initialSetUp(meins, prepro, boundaries, xmax, ymax, params,numOfCPUs);
    initialSetUp(meins, prepro, bound, xmax, ymax, params, numOfCPUs);
    write_vtk_output2D(meins, 0);
        
    
    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    //std::chrono::high_resolution_clock::duration parallel, sequential;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();

    const int outputInterval = timetrack.getOutputInt();

    int i;

    while (timetrack.proceed() == true)
    {
        start = std::chrono::high_resolution_clock::now();

        meins.collideAll(numOfCPUs,false,true);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        end = std::chrono::high_resolution_clock::now();
        parallel +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%outputInterval == 0) 
        {
            cout << i <<endl;
            write_vtk_output2D(meins, i);
        } 
    }

      cout<<"\nBerechnung beendet nach "<< sequential+parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon parallel:  "<< parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon sequentiell "<< sequential <<" Micro-Sekunden"<<endl;

    return 0;
}

void initialSetUp(Lattice2D& meins, Preprocess_Drop prepro, Boundaries bound, int xmax, int ymax, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();
    
    const Cell2D drop1(eqDistro({{rho_liquid,0}}, {{Vector2D(prepro.getUD(),0),Vector2D(0,0)}}, params.getPhi2D()));
    const Cell2D drop2(eqDistro({{rho_liquid,0}}, {{Vector2D(- prepro.getUD(),0),Vector2D(0,0)}}, params.getPhi2D()));
    const Cell2D gas(0,rho_gas,false);

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    // const int xm1 = xmax/2;
    const int xm1 = (xmax * 0.5) - 3 * R1;
    const int ym1 = ymax * 0.5 ;

    const int xm2 = (xmax * 0.5) + 3 * R1;
    const int ym2 = ymax * 0.5 ;
    //const int ym1 = R1 + xm1;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++)
        {
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,drop1);
            else if ( (i-xm2)*(i-xm2) + (j-ym2)*(j-ym2) < R1*R1 ) meins.setCell(i,j,drop2); 
            else meins.setCell(i,j,gas);
        }
    }

    meins.setBoundaries(bound);
  
    meins.equilibriumIni();

    cout<<"Initialisierung beendet\n"<<endl;
}

const Boundaries getBoundaries(double gamma)
{
    const BoundaryInformation generic = BoundaryInformation(2, ColSet({{0,1.0/gamma}}), VeloSet3D());

    Boundaries bound;
    bound.north = generic;
    bound.south = generic;
    bound.west = generic;
    bound.east = generic;

    return bound;
}