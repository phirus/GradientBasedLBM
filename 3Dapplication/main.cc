#include<iostream>
#include<ctime>
#include<vector>
#include<chrono>

#include"../src/3D/Lattice3D.h"
#include"../src/3D/BinaryIO3D.h"
#include"../src/3D/Analyze3D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice3D& meins, Preprocess& prepro, Boundaries& bound, ParamSet params, int numOfCPUs);
void initializeShearfFlow(Lattice3D& meins, Preprocess& prepro, Boundaries& bound, ParamSet params, int numOfCPUs);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("boundary,b", boost::program_options::value<string> (), "specify boundary input file")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ("override,o", boost::program_options::value<string> (), "specify parameter file to bypass preprocess routine")
        ("restart,r", boost::program_options::value<string> (), "specify restart file")
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
    Lattice3D bubble_only(xmax,ymax,zmax);

    if (vm.count("restart")) 
    {
        cout << "Restart file is: " << vm["restart"].as<string>() << ".\n" << endl ;
        Lattice3D tmpL;
        Preprocess tmpP;
        Timetrack tmpT;
        bool tmpB = read_restart_file3D(tmpL, tmpP, tmpT, vm["restart"].as<string>());
        if (tmpB == true)
        {
            meins = tmpL;
            prepro = tmpP;
            timetrack = tmpT;
        }
        else 
        {
            cout << "failed to read input file" << endl;
            return 1;
        }
        if (vm.count("preprocess"))
        {
            timetrack.setMaxCount( timetrack_input.getMaxCount() );
            timetrack.setOutputInt( timetrack_input.getOutputInt() );
            timetrack.setRestartInt( timetrack_input.getRestartInt() );
        }
    }
    else {      // vm.count("restart")
        initialSetUp(bubble_only, prepro, boundaries, params, numOfCPUs);
        initializeShearfFlow(meins, prepro, boundaries, params, numOfCPUs);

        meins.copyCellsFromOther(bubble_only, bubble_only.findBubbleCells());
        write_vtk_output3D(meins, 0);        
    }

    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();


    const int outputInterval = timetrack.getOutputInt();
    const int restartInterval = timetrack.getRestartInt();

    write_file_header("BubblePlot.csv", "time;PosX;PosY;PosZ;v_x;v_y;v_z");
    //write_file_header("BubblePosPlot.dat", "time \t PosX \t PosY \t PosZ");
    //write_file_header("BubbleVeloPlot.dat", "time \t vx \t vy \t vz");

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
        
        if(i%100 == 0) 
        {
            cout << i <<endl;
        }
        
        if(i%10 == 0) 
        {
            const boost::array<Vector3D,2> bubble_data = meins.getBubbleData();
            const Vector3D pos = bubble_data[0];
            const Vector3D velo = bubble_data[1];

            // const double Re = getReynolds(params, velo.z, prepro.getResolution());

            BubbleBox3D bubblebox = meins.getBubbleBox();
            bubblebox.setBubble(pos.x,pos.y,pos.z);
            meins.setBubbleBox(bubblebox);

            writeBubbleFitData(meins, createFilename("bubbleFit_", i, ".csv"));            
            //write_data_plot_linewise(i ,pos.x, pos.y , "BubblePosPlot.dat");
            //write_data_plot_linewise(i ,velo.x, velo.y, "BubbleVeloPlot.dat");
            write_csv_linewise(i, pos.x, pos.y, pos.z, velo.x, velo.y, velo.z, "BubblePlot.csv");
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output3D(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file3D(meins, prepro, timetrack, restart_file_name);
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

   for (int i = 1; i< 101; i++){
       meins.collideAll(numOfCPUs,false,false);
       meins.evaluateBoundaries();
       meins.streamAll(numOfCPUs);
       if(i%100 == 0) cout << i<<endl;
   }

   BubbleBox3D bubblebox;
   bubblebox.setBubble(xm,ym,zm);
   bubblebox.setDX(3*prepro.getResolution());
   bubblebox.setDY(3*prepro.getResolution());
   bubblebox.setDZ(3*prepro.getResolution());
   meins.setBubbleBox(bubblebox);


    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}

void initializeShearfFlow(Lattice3D& meins, Preprocess& prepro,  Boundaries& bound, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    const int xmax = prepro.getXCells();
    const int ymax = prepro.getYCells();
    const int zmax = prepro.getZCells();

    // get densities
    const double rho_liquid = prepro.getRhoL();

    const Cell3D liquid(rho_liquid,0,false);

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++)
        {
            for(int k=0;k<zmax;k++)
            {
                meins.setCell(i,j,k,liquid);
            }
        }
    }

    double grad = (bound.east.getVelocity()[0].z - bound.west.getVelocity()[0].z) / double(ymax);
    meins.setShearProfile(grad , bound.west.getVelocity()[0].z);
    meins.setBoundaries(bound);
    //meins.equilibriumIni();

    for (int i = 1; i< 101; i++)
    {
        meins.collideAll(numOfCPUs,false,false);
        meins.evaluateBoundaries();
        meins.streamAll(numOfCPUs);
   }


    cout<<"Shear flow steadily initialized\n"<<endl;
}