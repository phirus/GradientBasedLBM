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


    // create a Lattice   
    int ymax = prepro.getYCells();
    int xmax = prepro.getXCells();
    
    Lattice2D meins(xmax,ymax);
    initialSetUp(meins, prepro, boundaries, xmax, ymax, params, numOfCPUs);
    
    write_vtk_output2D(meins, 0);
    
    const int outputInterval = timetrack.getOutputInt();

    write_file_header("BubblePosPlot.dat", "time \t PosX \t PosY");

    bool isAppendOkay = true;


    int i;

    while (timetrack.proceed() == true)
    {
        meins.collideAll(numOfCPUs,true,true);
        meins.evaluateBoundaries(numOfCPUs);
        meins.streamAll(numOfCPUs);

        timetrack.timestep();

        // Output if necessary
        i = timetrack.getCount();
        cout << i << endl;
        
        if(i%outputInterval == 0) 
        {
            write_vtk_output2D(meins, i);
        } 
        
       
        if(i%5 == 0) 
        {
            const boost::array<Vector2D,3> bubble_data = meins.getBubbleData();
            const Vector2D pos = bubble_data[0];

            BubbleBox2D bubblebox = meins.getBubbleBox();
            bubblebox.setBubble(pos.x,pos.y);
            meins.setBubbleBox(bubblebox);

            write_data_plot_linewise(i ,pos.x, pos.y + meins.getOffset(), "BubblePosPlot.dat");

            // if (pos.y > 0.5 * ymax && isCutoffOkay == true)
            // {
            //     isCutoffOkay = false;
            //     meins = meins.latticeCutOff(ymax / 10);
            // }

            if (pos.y > 600 && isAppendOkay == true)
            {
//                write_vtk_output2D(meins, i);
                meins = meins.latticeAppend(1000,prepro.getRhoL(),0);
  //              write_vtk_output2D(meins, i+1);
                isAppendOkay = false;
            }

        }         
    }

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
    const int xm1 = xmax * 0.5;
    const int ym1 = 2*R1;
    
    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            // if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // // if( j > ym1) meins.setCell(i,j,air);
            // else meins.setCell(i,j,liquid);
            meins.setCell(i,j,liquid);
        }
    }

    meins.setBoundaries(bound);
  
    meins.equilibriumIni();

   for (int i = 1; i< 11; i++)
   {
       meins.collideAll(numOfCPUs,false,false);
       meins.evaluateBoundaries();
       meins.streamAll(numOfCPUs);
   }

    BubbleBox2D bubblebox;
    bubblebox.setBubble(xm1,ym1);
    bubblebox.setW(5*prepro.getResolution());
    bubblebox.setH(5*prepro.getResolution());
    meins.setBubbleBox(bubblebox);

    cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
}


