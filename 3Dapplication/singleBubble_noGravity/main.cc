#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../../src/3D/Lattice3D.h"
#include"../../src/3D/BinaryIO3D.h"
#include"../../src/3D/Analyze3D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice3D& meins, Preprocess& prepro, Boundaries& bound, ParamSet params, int numOfCPUs);

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
    int xmax = prepro.getXCells();
    int ymax = prepro.getYCells();
    int zmax = prepro.getZCells();

    Lattice3D meins(xmax,ymax,zmax);

    initialSetUp(meins, prepro, boundaries, params, numOfCPUs);
    write_vtk_output3D_verboose(meins, 0);        

    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();


    const int outputInterval = timetrack.getOutputInt();

    write_file_header("BubblePlot.csv", "time;PosX;PosY;PosZ;v_x;v_y;v_z");

    int i;

    while (timetrack.proceed() == true)
    {
        start = std::chrono::high_resolution_clock::now();

        meins.collideAll(numOfCPUs);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        end = std::chrono::high_resolution_clock::now();
        parallel +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

        start = std::chrono::high_resolution_clock::now();
        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%10 == 0) 
        {
            cout << i <<endl;
        }
        
        if(i%50 == 0) 
        {
            const boost::array<Vector3D,2> bubble_data = meins.getBubbleData();
            const Vector3D pos = bubble_data[0];
            const Vector3D velo = bubble_data[1];

            BubbleBox3D bubblebox = meins.getBubbleBox();
            bubblebox.setBubble(pos.x,pos.y,pos.z);
            meins.setBubbleBox(bubblebox);

            writeBubbleFitData(meins, createFilename("bubbleFit_", i, ".csv"));            
            write_csv_linewise(i, pos.x, pos.y, pos.z, velo.x, velo.y, velo.z, "BubblePlot.csv");
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output3D_verboose(meins, i);
        } 
        
        end = std::chrono::high_resolution_clock::now();
        sequential +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();         
    }

    cout<<"\nBerechnung beendet nach "<< sequential+parallel <<" Micro-Sekunden"<<endl;
    cout<<"\nDavon parallel:  "<< parallel <<" Micro-Sekunden"<<endl;
    cout<<"\nDavon sequentiell "<< sequential <<" Micro-Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice3D& meins, Preprocess& prepro, Boundaries& bound, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell3D air(0,rho_gas,false);
    const Cell3D liquid(rho_liquid,0,false);

    //const Cell3D wall(0,0,true);

    const int xmax = prepro.getXCells();
    const int ymax = prepro.getYCells();
    const int zmax = prepro.getZCells();

    // setup geometry (bubble at the bottom, x and y centered)
    const int R = prepro.getResolution()/2;
    const int xm = xmax/2;
    const int ym = ymax/2;
    const int zm = zmax/2;
    // const int zm =  R + 20;

    for(int i=0; i< xmax; i++)
    {
        for(int j=0; j< ymax; j++)
        {
            for(int k=0; k< zmax; k++)
            {
                if( (i-xm)*(i-xm) + (j-ym)*(j-ym) + (k-zm)*(k-zm) < R*R ) meins.setCell(i,j,k,air);
                // if( k > zm) meins.setCell(i,j,air);
                else meins.setCell(i,j,k,liquid);
            }
        }
    }

    meins.equilibriumIni();

   BubbleBox3D bubblebox;
   bubblebox.setBubble(xm,ym,zm);
   bubblebox.setDX(3*prepro.getResolution());
   bubblebox.setDY(3*prepro.getResolution());
   bubblebox.setDZ(3*prepro.getResolution());
   meins.setBubbleBox(bubblebox);

    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}
