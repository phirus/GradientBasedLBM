#include<iostream>
#include<ctime>
#include<vector>

#include"../src/Lattice.h"
#include"../src/BinaryIO.h"
#include"../src/Analyze.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax, ParamSet params);

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
    int ymax = prepro.getHeight();
    int xmax = prepro.getWidth();
    Lattice meins(xmax,ymax);

    if (vm.count("restart")) 
    {
        cout << "Restart file is: " << vm["restart"].as<string>() << ".\n" << endl ;
        Lattice tmpL;
        Preprocess tmpP;
        Timetrack tmpT;
        bool tmpB = read_restart_file(tmpL, tmpP, tmpT, vm["restart"].as<string>());
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
        initialSetUp(meins, prepro, xmax, ymax, params);
        write_vtk_output(meins, 0);;
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
            write_vtk_output(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file(meins, prepro, timetrack, restart_file_name);
        }         
    }

    time(&end);
    write_data_plot(reynolds_data, 1000, "ReynoldsPlot.dat");
    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell air(0,rho_gas,false);
    const Cell liquid(rho_liquid,0,false);

    const Cell wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax/2;
    // const int ym1 = 2*R1;
    const int ym1 = R1 + 20;

    // const int ym1 = ymax/2;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.bottomWall();
    meins.equilibriumIni();

    // write_techplot_output(meins,0,true);
    // // temporal mass_balance
    // std::vector<double> count;
    // std::vector<double> liquid_mass;
    // std::vector<double> gas_mass;

   for (int i = 1; i< 1001; i++){
       meins.collideAll(1,false,false);
       meins.streamAll(1);
       if(i%100 == 0) cout << i<<endl;
       // if(i%10 == 0) write_techplot_output(meins,i,true);;
   }

   // for (int i = 501; i< 100001; i++){
    
   //     meins.collideAll(4,false,false);
   //     meins.streamAll(4);
   //     if(i%500 == 0)
   //     {
   //      double tmpL,tmpG;
   //      meins.mass_balance(tmpL,tmpG);
   //      count.push_back(i);
   //      liquid_mass.push_back(tmpL);
   //      gas_mass.push_back(tmpG);
   //     }
   // }
   //        meins.collideAll(4,false,false,true);
    // write_data_plot(count, liquid_mass, gas_mass);
    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}