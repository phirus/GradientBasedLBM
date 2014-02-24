#include<iostream>
#include<ctime>
#include<vector>

#include"../src/Lattice.h"
#include"../src/BinaryIO.h"
#include<boost/program_options.hpp>

using namespace std;
namespace po = boost::program_options;

void explicit_output(Lattice& l, int x, int y);

void plot_array(const array& a);
void plot_distribution(const DistributionSetType& f);
void plot_cell(const Cell& c);

int main(int argc, char** argv){

    boost::program_options::options_description desc("Allowed options");
	desc.add_options()
        ("help,h", "produce help message")
        ("preprocess,p", boost::program_options::value<string> (), "specify preprocess parameter file")
        // ("restart,r", boost::program_options::value<string> (), "specify restart file")
        ;
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc),vm);
    boost::program_options::notify(vm);

    if(vm.count("help")){
        cout << desc << endl;
        return 1;
    }

    Preprocess prepro = read_preprocess_file("preprocessFile");
    Timetrack timetrack = read_timetrack_file(prepro, "preprocessFile");

    if (vm.count("preprocess")) {
        cout << "preprocess file is: " << vm["preprocess"].as<string>() << ".\n" << endl ;
        prepro = read_preprocess_file(vm["preprocess"].as<string>());
        timetrack = read_timetrack_file(prepro, vm["preprocess"].as<string>());
    }
   
    int ymax = 100;
    int xmax = 100;

    Lattice meins(xmax,ymax);
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
    const int ym1 = ymax/2;

    for(int j=0; j< ymax; j++)
    {
        for(int i=0; i< xmax; i++){
            if( (i-xm1)*(i-xm1) + (j-ym1)*(j-ym1) < R1*R1 ) meins.setCell(i,j,air);
            // if( j > ym1) meins.setCell(i,j,air);
            else meins.setCell(i,j,liquid);
        }
    }

    meins.equilibriumIni();
    write_techplot_output(meins,0,true);

    // Bildung der Grenzschicht bevor Schwerkraft zugeschaltet wird

   // for (int i = 1; i< 1001; i++){
    
   //     meins.collideAll(4,false,false);
   //     meins.streamAll(4);

    
   //     if(i%100 == 0) cout << i<<endl;
   //     if(i%500 == 0) write_techplot_output(meins,i,true);
   // }
    // meins.collideAll(4,false,false,true);
   explicit_output(meins,1,1);

   array a;
   plot_array(a);


   return 0;
}

void explicit_output(Lattice& l, int x, int y)
{
    cout<<"Test"<<endl;
}

void plot_array(const array& a)
{
    for(int i=0; i<9;i++){
        cout<<"\n"<< a[i];
    }
}

void plot_distribution(const DistributionSetType& f)
{
    cout << "\ndense phase:";
    plot_array(f[0]);
    cout << "\n\ndilute phase:";
    plot_array(f[1]);
}

void plot_cell(const Cell& c)
{
    const ColSet rho = c.getRho();
    cout << "\nrho[0] = " << rho[0];
    cout << "\nrho[1] = " << rho[1];
    cout << "\ndeltaRho = " << c.getDeltaRho();

    const VeloSet u = c.getU();
    cout << "\nu[0]: x = " << u[0].x << "\ty=" << u[0].y;
    cout << "\nu[1]: x = " << u[1].x << "\ty=" << u[1].y;       

    cout << "\n";

    plot_distribution(c.getF());
}