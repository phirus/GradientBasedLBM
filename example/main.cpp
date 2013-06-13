#include <iostream>
#include<ctime>

#include"../src/lattice.h"
#include"../src/binaryIO.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

int main(int argc, char** argv){

    time_t start,end;
	time(&start);

	boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);


    if(vm.count("help")){
        cout << desc << endl;
        return 1;
    }

    int ymax = 150;
    int xmax = 80;

//    double height = 1; // /m
//    double dx = height/ymax;
    double dx = 5e-5;

    double rho_liquid = 1;
    double rho_gas= 0.001;
    double ratio = rho_liquid / rho_gas;

    Cell air(0,rho_gas,false);
    Cell liquid(rho_liquid,0,false);
    Cell wall(0,0,true);

    ParamSet params;
    params.setAlpha(0.2);
    params.setRatio(rho_liquid,ratio);
    params.setSigma(6.74e-15);
    params.setBeta(0.9);

    params.setDeltaX(dx);
    params.setSoundSpeed(3);
//    params.setOmega(1.9597,1.9597,0.1);

    Lattice meins(xmax,ymax);
    meins.setParams(params);

    int R1 = 20;

    int xm1 = xmax/2;
    int ym1 = 2*R1;//-R1;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }


//    int R2 = R1;
//
//    int xm2 = xmax/2;
//    int ym2 = ymax/2+R1;
//
//    for(int j=0; j< ymax; j++)
//    {
//        for(int i=0; i< xmax; i++){
//            if( (i-xm2)*(i-xm2) + (j-ym2)*(j-ym2) < R2*R2 ) meins.setCell(i,j,air);
//        }
//    }

//    int height = 40;
//    int width = 40;
//
//    for(int j=0; j< ymax; j++)
//    {
//        for(int i=0; i< xmax; i++){
//            if( i > (xm - width/2) && i < (xm + width/2)   && j > (ym - height/2) && j < (xm + height/2)  ) meins.setCell(j,i,liquid);
//            else meins.setCell(j,i,air);
//        }
//    }


    meins.bottomWall();
    meins.equilibriumIni();


//    for (int i = 1; i< 501; i++){
//        meins.collideAll(1,false);
//        meins.streamAll(1);
//        if(i%50 == 0) cout << i<<endl;
////        if(i%1000 == 0)  meins.techplotOutput(i,true);
//    }
//
//cout<<"Initialisierung beendet\n\nSchwerkraft wird zugeschaltet\n"<<endl;
////


//int numSteps = 1e6 + 1;
//for (int i = 0; i< numSteps; i++){
//        meins.collideAll(4);
//        meins.streamAll(4);
//        if(i%1000 == 0) cout << i<<endl;
//        if(i%10000 == 0)  techplotOutput(meins,i,true);
//        if(i%10000 == 0) binary_output(meins);
//    }




    time(&end);
    cout<<"\nBerechnung beendet nach "<< difftime(end,start) <<" Sekunden"<<endl;

    return 0;
}
