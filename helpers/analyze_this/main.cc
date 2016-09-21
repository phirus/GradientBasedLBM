#include<iostream>
#include<ctime>
#include<vector>

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

    Preprocess prepro = read_preprocess_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file("preprocessFile");
    ParamSet params = prepro.getParamSet();
    Boundaries boundaries = read_boundaries_file("BoundaryInput");

    if (vm.count("override")) 
    {
        cout << "override preprocess file with: " << vm["override"].as<string>() << ".\n" << endl ;
        ParamSet params_input = read_paramset_file(vm["override"].as<string>());
        params = params_input;
    }

    // create a Lattice   
    int ymax = prepro.getYCells();
    int xmax = prepro.getXCells();
    Lattice2D meins(xmax,ymax);
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();
    const double rho_gas = rho_liquid / prepro.getGamma();

    const Cell2D air(0,rho_gas,false);
    const Cell2D liquid(rho_liquid,0,false);


    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax * 0.5;
    const int ym1 = R1 + 20;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.setBoundaries(boundaries);
    meins.equilibriumIni();

   for (int i = 1; i< 101; i++){
       meins.collideAll(1,false,false);
       //meins.evaluateBoundaries();
       meins.streamAll(1);
   }

    write_vtk_output2D(meins, 0);    

    time_t start,end;
    time(&start);

    const int outputInterval = timetrack.getOutputInt();
    const int restartInterval = timetrack.getRestartInt();

    std::vector<double> iter_count_data;
    std::vector<double> bubble_pos_x_data;
    std::vector<double> bubble_pos_y_data;
    std::vector<double> x_velo_data;
    std::vector<double> y_velo_data;
    std::vector<double> reynolds_data;
    
    int i;

    while (timetrack.proceed() == true)
    {
        meins.collideAll(1,true,true);
        meins.evaluateBoundaries();
        meins.streamAll(1);

        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        cout << i <<endl;

        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file2D(meins, prepro, timetrack, restart_file_name);
        }
        
        iter_count_data.push_back(i);

        const double reynolds_tmp = getReynolds(meins, prepro.getResolution());            
        reynolds_data.push_back(reynolds_tmp);
        write_data_plot(reynolds_data, 1, "ReynoldsPlot.dat");

        const Vector2D velo_tmp = getBubbleVelocity(meins);            
        x_velo_data.push_back(velo_tmp.x);
        y_velo_data.push_back(velo_tmp.y);

        write_data_plot(y_velo_data, x_velo_data, 1, "BubbleVeloPlot.dat");

        const Vector2D pos_tmp = getBubblePosition(meins);

        bubble_pos_x_data.push_back(pos_tmp.x);
        bubble_pos_y_data.push_back(pos_tmp.y);

        write_data_plot(bubble_pos_x_data, bubble_pos_y_data, 1, "BubblePosPlot.dat");

        nested_vector statistics;
        statistics.push_back(iter_count_data);
        statistics.push_back(bubble_pos_x_data);
        statistics.push_back(bubble_pos_y_data);
        statistics.push_back(x_velo_data);
        statistics.push_back(y_velo_data);
        statistics.push_back(reynolds_data);
        write_csv(statistics, "BubblePlot.csv", "time;PosX;PosY;v_x;v_y;Re");
    }

    time(&end);
    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;

    return 0;
}