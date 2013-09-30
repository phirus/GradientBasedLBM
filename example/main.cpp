#include <iostream>
#include<ctime>

#include"../src/lattice.h"
#include"../src/binaryIO.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("cpu,c", boost::program_options::value<int> (), "takes the number of CPUs")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        ("restart,r", boost::program_options::value<string> (), "specify restart file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help")){
        cout << desc << endl;
        return 1;
    }

    int numOfCPUs = 1;
    Preprocess prepro = getFilePreprocess("preprocessFile");
    Timetrack timetrack = getFileTimetrack(prepro, "preprocessFile");    

    if (vm.count("preprocess")) {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro = getFilePreprocess(vm["preprocess"].as<string>());
        Timetrack timetrack = getFileTimetrack(prepro, vm["preprocess"].as<string>());
    }
   
    int ymax = 150;
    int xmax = 80;
    // create a Lattice
    Lattice meins(xmax,ymax);
    initialSetUp(meins, prepro, xmax, ymax);
       
    if (vm.count("restart")) {
        cout << "Restart file is: " << vm["restart"].as<string>() << ".\n" << endl ;
        Lattice tmpL;
        Preprocess tmpP;
        Timetrack tmpT;
        bool tmpB = restart_read(tmpL, tmpP, tmpT, vm["restart"].as<string>());
        if (tmpB == true){
            meins = tmpL;
            prepro = tmpP;
            timetrack = tmpT;
        }
        else {
            cout << "failed to read input file" << endl;
            return 1;
        }
    }

    if(vm.count("cpu")){
        numOfCPUs = vm["cpu"].as<int>();
        cout << "number of CPUs set to "<< numOfCPUs;
    }

    time_t start,end;
    time(&start);

    while (timetrack.proceed() == true){
        try{
            meins.collideAll(numOfCPUs);
            meins.streamAll(numOfCPUs);
            timetrack.timestep();
            int i = timetrack.getCount();
            if(i%1000 == 0) cout << i<<endl;
            if(i%10000 == 0)  techplotOutput(meins,i,true);
            if(i%10000 == 0) restart_file(meins, prepro, timetrack);
        }
        catch(string s)
        {
            prepro.refine();
            const ParamSet params = prepro.getParamSet();
            meins.setParams(params);
            timetrack.refine();
            cout << s << endl;
            cout<<"\nGitter verfeinert bei i = " << timetrack.getCount() << endl;
        }
        
    }

    time(&end);
    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;

    return 0;
}


void initialSetUp(Lattice& meins, Preprocess& prepro, int xmax, int ymax){
    // set the parameters    
    const ParamSet params = prepro.getParamSet();
    meins.setParams(params);

    // get densities
    const double rho_liquid = prepro.convertRhoL();
    const double rho_gas = prepro.convertRhoG();

    const Cell air(0,rho_gas,false);
    const Cell liquid(rho_liquid,0,false);
    const Cell wall(0,0,true);

    // setup geometry (bubble at the bottom, x-centered)
    const int R1 = prepro.getResolution()/2;
    const int xm1 = xmax/2;
    const int ym1 = 2*R1;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.bottomWall();
    meins.equilibriumIni();

    // Bildung der Grenzschicht bevor Schwerkraft zugeschaltet wird

   for (int i = 1; i< 501; i++){
       meins.collideAll(1,false);
       meins.streamAll(1);
       if(i%100 == 0) cout << i<<endl;
   }
cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
//

}