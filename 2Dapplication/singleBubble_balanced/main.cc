#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../../src/2D/Lattice2D.h"
#include"../../src/2D/BinaryIO2D.h"
#include"../../src/2D/Analyze2D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice2D& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params, int numOfCPUs);
const Boundaries setVelocity(const Boundaries& bound, const Vector3D& v_inc);


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
    //int ymax = prepro.getYCells();
    int ymax = prepro.getXCells();
    int xmax = prepro.getXCells();
    Lattice2D meins(xmax,ymax);

    initialSetUp(meins, prepro, boundaries, xmax, ymax, params,numOfCPUs);
    write_vtk_output2D(meins, 0);
    

    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();

    const int outputInterval = timetrack.getOutputInt();

    write_file_header("BubblePlot.csv", "time;PosX;PosY;v_x;v_y;");
    write_file_header("Massbalance.csv", "time;liquidMass;gasMass;");
    write_file_header("BoundaryPlot.csv", "time;v_x;v_y;");
    int i;

    while (timetrack.proceed() == true)
    {
        start = std::chrono::high_resolution_clock::now();

        meins.collideAll(numOfCPUs,true);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        end = std::chrono::high_resolution_clock::now();
        parallel +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

        start = std::chrono::high_resolution_clock::now();
        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%100 == 0) 
        {
            cout << i <<endl;
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        } 
       
        if(i%10 == 0) 
        {
            const boost::array<Vector2D,3> bubble_data = meins.getBubbleData();
            const Vector2D pos = bubble_data[0];
            const Vector2D velo = bubble_data[1];

            const Boundaries bound_tmp = setVelocity(meins.getBoundaries(), Vector3D(0,velo.y,0));
            meins.setBoundaries(bound_tmp);
            write_data_plot_linewise(i , bound_tmp.north.getVelocity()[0].x, bound_tmp.north.getVelocity()[0].y,"BoundaryPlot.csv");


            double liquid_mass, gas_mass;
            meins.mass_balance(liquid_mass, gas_mass);
            write_data_plot_linewise(i , liquid_mass, gas_mass,"Massbalance.csv");           
            
            BubbleBox2D bubblebox = meins.getBubbleBox();
            bubblebox.setBubble(pos.x,pos.y);
            meins.setBubbleBox(bubblebox);

            writeBubbleFitData(meins, createFilename("bubbleFit_", i, ".csv"));            
            write_csv_linewise2D(i, pos.x, pos.y, velo.x, velo.y, "BubblePlot.csv");
        }

        end = std::chrono::high_resolution_clock::now();        
        sequential +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count(); 
    }
      cout<<"\nBerechnung beendet nach "<< sequential+parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon parallel:  "<< parallel <<" Micro-Sekunden"<<endl;
      cout<<"\nDavon sequentiell "<< sequential <<" Micro-Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice2D& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell2D air(0,rho_gas,false);
    const Cell2D liquid(rho_liquid,0,false);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax/2;
    const int ym1 = ymax/2;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.equilibriumIni();

    for (int i = 1; i< 101; i++)
    {
        meins.collideAll(numOfCPUs,false);
        meins.evaluateBoundaries();
        meins.streamAll(numOfCPUs);
        //if(i%100 == 0) cout << i<<endl;
    }

    meins.setBoundaries(bound);
  
    BubbleBox2D bubblebox;
    bubblebox.setBubble(xm1,ym1);
    bubblebox.setW(3*prepro.getResolution());
    bubblebox.setH(3*prepro.getResolution());
    meins.setBubbleBox(bubblebox);

    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}

const Boundaries setVelocity(const Boundaries& bound, const Vector3D& v_inc)
{
    Vector3D v = bound.north.getVelocity()[0] - (v_inc * 0.05) ;
    Vector3D zero(0,0,0);

    Boundaries b = bound;

    b.north.setVelocity({{v,zero}});
    b.east.setVelocity({{v,zero}});
    b.west.setVelocity({{v,zero}});
    b.south.setVelocity({{v,zero}});

    return b;
}