#include<iostream>
#include<ctime>
#include<vector>

#include"../src/3D/Lattice3D.h"
#include"../src/3D/BinaryIO3D.h"
#include"../src/3D/Analyze3D.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice3D& meins, Preprocess& prepro, int xmax, int ymax, int zmax, ParamSet params, int numOfCPUs);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ("bypass,b", boost::program_options::value<string> (), "specify parameter file to bypass preprocess routine")
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

    if (vm.count("bypass")) 
    {
        cout << "bypass preprocess file with: " << vm["bypass"].as<string>() << ".\n" << endl ;
        params_input = read_paramset_file(vm["bypass"].as<string>());
        params = params_input;
    }

    // create a Lattice   
    int xmax = prepro.getXCells();
    int ymax = prepro.getYCells();
    int zmax = prepro.getZCells();
    Lattice3D meins(xmax,ymax,zmax);

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
        initialSetUp(meins, prepro, xmax, ymax, zmax, params, numOfCPUs);
        write_vtk_output3D(meins, 0);;
    }

    time_t start,end;
    time(&start);

    const int outputInterval = timetrack.getOutputInt();
    const int restartInterval = timetrack.getRestartInt();

    std::vector<double> reynolds_data;
    int i;

    while (timetrack.proceed() == true)
    {
        meins.collideAll(numOfCPUs,true,true);
        meins.streamAll(numOfCPUs);
        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%10000 == 0) 
        {
            cout << i <<endl;
        }
        
        if(i%10 == 0) 
        {
            const double reynolds_tmp = getReynolds(meins, prepro.getResolution());            
            // reynolds_data.push_back(getReynolds(meins, prepro.getResolution()));
            reynolds_data.push_back(reynolds_tmp);
            write_data_plot(reynolds_data, 10, "ReynoldsPlot.dat");
            if(reynolds_tmp < 0) 
            {
                cout <<"\nReynolds < 0, probably reached the top "<<endl;
                break;
            }
        }

        if(i%outputInterval == 0) 
        {
            // write_techplot_output(meins,i,true);
            write_vtk_output3D(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file3D(meins, prepro, timetrack, restart_file_name);
        }         
    }

    time(&end);
    write_data_plot(reynolds_data, 1000, "ReynoldsPlot.dat");
    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice3D& meins, Preprocess& prepro, int xmax, int ymax, int zmax, ParamSet params, int numOfCPUs)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell3D air(0,rho_gas,false);
    const Cell3D liquid(rho_liquid,0,false);

    const Cell3D wall(0,0,true);

    // setup geometry (bubble at the bottom, x and y centered)
    const int radius = prepro.getResolution()/2;
    const int xm = xmax/2;
    const int ym = ymax/2;
    const int zm =  radius + 20;

    for(int i=0; i< xmax; i++)
    {
        for(int j=0; j< ymax; j++)
        {
            for(int k=0; k< zmax; k++)
            {
                if( (i-xm)*(i-xm) + (j-ym)*(j-ym) + (k-zm)*(k-zm) < radius*radius ) meins.setCell(i,j,k,air);
                // if( k > zm) meins.setCell(i,j,air);
                else meins.setCell(i,j,k,liquid);
            }
        }
    }

    meins.bottomWall();
    meins.equilibriumIni();

   for (int i = 1; i< 1001; i++){
       meins.collideAll(numOfCPUs,false,false);
       meins.streamAll(numOfCPUs);
       if(i%100 == 0) cout << i<<endl;
   }

    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}