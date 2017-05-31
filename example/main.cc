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

void initialSetUp(Lattice2D& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params, int numOfCPUs);
void initializeShearfFlow(Lattice2D& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params);

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
    //int ymax = prepro.getYCells();
    int ymax = 1000;
    int xmax = prepro.getXCells();
    Lattice2D meins(xmax,ymax);
    Lattice2D bubble_only(xmax,ymax);

    if (vm.count("restart")) 
    {
        cout << "Restart file is: " << vm["restart"].as<string>() << ".\n" << endl ;
        Lattice2D tmpL;
        Preprocess tmpP;
        Timetrack tmpT;
        bool tmpB = read_restart_file2D(tmpL, tmpP, tmpT, vm["restart"].as<string>());
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
    else {      // vm.count("restart"), if no restart file -> initialize

        //initialSetUp(meins, prepro, boundaries, xmax, ymax, params,numOfCPUs);
        
        initialSetUp(bubble_only, prepro, boundaries, xmax, ymax, params,numOfCPUs);
        initializeShearfFlow(meins, prepro, boundaries, xmax, ymax, params);
        write_vtk_output2D(meins, 0);
        
        meins.copyCellsFromOther(bubble_only, bubble_only.findBubbleCells());
        meins.setBubbleBox(bubble_only.getBubbleBox());
        write_vtk_output2D(meins, 0);
    }

    //time_t start,end;
    std::chrono::high_resolution_clock::time_point start,end;
    //std::chrono::high_resolution_clock::duration parallel, sequential;
    auto parallel = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();
    auto sequential = std::chrono::duration_cast<std::chrono::microseconds>( end - end).count();

    const int outputInterval = timetrack.getOutputInt();
    const int restartInterval = timetrack.getRestartInt();

    write_file_header("BubblePlot.csv", "time;PosX;PosY;v_x;v_y;Re;Fx;Fy");
    write_file_header("BubblePosPlot.dat", "time \t PosX \t PosY");
    write_file_header("BubbleVeloPlot.dat", "time \t vx \t vy");
    write_file_header("ForcePlot.dat", "time \t Fx \t Fy");

    int i;

//    bool isCutoffOkay = true;
    bool isAppendOkay = true;

    while (timetrack.proceed() == true)
    {
        start = std::chrono::high_resolution_clock::now();

        meins.collideAll(numOfCPUs,true,true);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        end = std::chrono::high_resolution_clock::now();
        parallel +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

        start = std::chrono::high_resolution_clock::now();
        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        
        if(i%10000 == 0) 
        {
            cout << i <<endl;
        }

        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        } 
        
        if(i%restartInterval == 0)
        {
            const string restart_file_name =  createFilename("restart", i, ".bin");
            write_restart_file2D(meins, prepro, timetrack, restart_file_name);
        }
        
        if(i%10 == 0) 
        {
            const boost::array<Vector2D,3> bubble_data = meins.getBubbleData();
            const Vector2D pos = bubble_data[0];
            const Vector2D velo = bubble_data[1];           
            const Vector2D force = bubble_data[2];

            const double Re = getReynolds(params, velo.y, prepro.getResolution());

            BubbleBox2D bubblebox = meins.getBubbleBox();
            bubblebox.setBubble(pos.x,pos.y);
            meins.setBubbleBox(bubblebox);

            writeBubbleFitData(meins, createFilename("bubbleFit_", i, ".csv"));            
            write_data_plot_linewise(i ,pos.x, pos.y + meins.getOffset(), "BubblePosPlot.dat");
            write_data_plot_linewise(i ,velo.x, velo.y, "BubbleVeloPlot.dat");
            write_csv_linewise(i, pos.x, pos.y + meins.getOffset(), velo.x, velo.y, Re, force.x, force.y, "BubblePlot.csv");

            write_data_plot_linewise(i ,force.x, force.y, "ForcePlot.dat");

            end = std::chrono::high_resolution_clock::now();
            sequential +=  std::chrono::duration_cast<std::chrono::microseconds>( end - start).count();

            // if (pos.y > 0.5 * ymax && isCutoffOkay == true)
            // {
            //     isCutoffOkay = false;
            //     meins = meins.latticeCutOff(ymax / 10);
            // }

            if (pos.y > 600 && isAppendOkay == true)
            {
                meins = meins.latticeAppend(1000,prepro.getRhoL(),0);
                isAppendOkay = false;
            }

            if(pos.y  > 0.95 * 2000)
            {
                cout << "\nBubble reached the top";
                break;
            }

            // if(pos.y + meins.getOffset() > 0.95 * ymax)
            // {
            //     cout << "\nBubble reached the top";
            //     break;
            // }
        }         
    }

//    time(&end);
//    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;
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

    // const Cell2D wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    // const int xm1 = xmax/2;
    const int xm1 = xmax * 0.5;
    const int ym1 = 2*R1;
    //const int ym1 = R1 + xm1;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    //meins.setBoundaries(bound);
  
    meins.equilibriumIni();

//   for (int i = 1; i< 1001; i++)
   for (int i = 1; i< 101; i++)
   {
       meins.collideAll(numOfCPUs,false,false);
       //meins.evaluateBoundaries();
       meins.streamAll(numOfCPUs);
       if(i%100 == 0) cout << i<<endl;
   }

    BubbleBox2D bubblebox;
    bubblebox.setBubble(xm1,ym1);
    bubblebox.setW(5*prepro.getResolution());
    bubblebox.setH(5*prepro.getResolution());
    meins.setBubbleBox(bubblebox);

    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}


void initializeShearfFlow(Lattice2D& meins, Preprocess& prepro, Boundaries& bound, int xmax, int ymax, ParamSet params)
{
    // set the parameters        
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.getRhoL();

    const Cell2D liquid(rho_liquid,0,false);

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            meins.setCell(i,j,liquid);
        }
    }

    meins.setBoundaries(bound);
    double grad = - bound.west.getVelocity()[0].y / double(xmax);
    meins.setShearProfile(grad , bound.west.getVelocity()[0].y);
    meins.equilibriumIni();

    cout<<"Shear flow steadily initialized\n"<<endl;
}
