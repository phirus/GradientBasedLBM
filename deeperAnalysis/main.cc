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
void plot_cell(const Cell& c, bool verbose);
void plot_gradient_information(Lattice& l, int x, int y);
const Cell mock_single_collision(const Cell& c, const ParamSet& params);

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
   for (int i = 1; i< 501; i++){
    
       meins.collideAll(1,false,false);
       meins.streamAll(1);
   }
   explicit_output(meins,37,37);

   return 0;
}

void explicit_output(Lattice& l, int x, int y)
{
    cout.precision(15);
    Cell tmp;
    tmp = l.getCell(x,y);
    cout << "\nactual cell:";
    plot_cell(tmp,true);

    cout << "\n";
    cout << "\n";

    const ParamSet param = l.getParams();
    const DistributionSetType phi = param.getPhi();
    const ColSet rho_k = tmp.getRho();
    const double rho = sum(rho_k);
    const VeloSet u = tmp.getU();
    const DistributionSetType fEq = eqDistro(rho_k, u, phi);

    Cell EqCell(fEq);
    EqCell.calcRho();
    cout << "\nequilibrium cell:";
    plot_cell(EqCell,true);

    cout << "\n";
    cout << "\n";

    cout << "\npost-collision cell:";
    cout << "\nomega = "<<param.getOmega(0);
    plot_cell(mock_single_collision(tmp, param),true);

    cout << "\n";
    cout << "\n";
}

void plot_array(const array& a)
{
    for(int i=0; i<9;i++){
        cout<<"\n"<< a[i] << "\tf[" << i << "]";
    }
}

void plot_distribution(const DistributionSetType& f)
{
    cout << "\ndense phase:";
    plot_array(f[0]);
    cout << "\n\ndilute phase:";
    plot_array(f[1]);
}

void plot_cell(const Cell& c, bool verbose)
{
    const ColSet rho = c.getRho();
    cout << "\nrho[0] = " << rho[0];
    cout << "\nrho[1] = " << rho[1];
    cout << "\ndeltaRho = " << c.getDeltaRho();
    if(true == verbose)
    {
        const VeloSet u = c.getU();
        cout << "\nu[0]: x = " << u[0].x << "\ty=" << u[0].y;
        cout << "\nu[1]: x = " << u[1].x << "\ty=" << u[1].y; 
    }
    cout << "\n";
    plot_distribution(c.getF());
}

void plot_gradient_information(Lattice& l, int x, int y)
{
    Vector gradient = l.getGradient(x,y);
    
    {
        cout.precision(15);
        cout << "\ndelta";
        cout << "\nn: "  << "\t" << l.getCell(x  ,y+1).getDeltaRho();
        cout << "\nne: " << "\t" << l.getCell(x-1,y+1).getDeltaRho();
        cout << "\ne: "  << "\t" << l.getCell(x-1,y  ).getDeltaRho();
        cout << "\nse: " << "\t" << l.getCell(x-1,y-1).getDeltaRho();
        cout << "\ns: "  << "\t" << l.getCell(x  ,y-1).getDeltaRho();
        cout << "\nsw: " << "\t" << l.getCell(x+1,y-1).getDeltaRho();
        cout << "\nw: "  << "\t" << l.getCell(x+1,y  ).getDeltaRho();
        cout << "\nnw: " << "\t" << l.getCell(x+1,y+1).getDeltaRho();
        cout << "\nnn: " << "\t" << l.getCell(x  ,y+2).getDeltaRho();
        cout << "\nee: " << "\t" << l.getCell(x-2,y  ).getDeltaRho();
        cout << "\nss: " << "\t" << l.getCell(x  ,y-2).getDeltaRho();
        cout << "\nww: " << "\t" << l.getCell(x+2,y  ).getDeltaRho();
    }
    
    cout << "\ngradient: x= " << gradient.x << "\ty = " << gradient.y ;
}

const Cell mock_single_collision(const Cell& c, const ParamSet& param)
{
    const DistributionSetType phi = param.getPhi();
    const RelaxationPar relax = param.getRelaxation();
    DistributionSetType  fTmp;
    const DistributionSetType fCell = c.getF();
    
    const ColSet rho_k = c.getRho();
    const VeloSet u = c.getU();
    const DistributionSetType fEq = eqDistro(rho_k, u, phi,true);
    const double omega = param.getOmega(c.calcPsi());
                
    for (int q=0; q<9; q++)
    {
        for (int color=0;color<=1; color++)
        {
            fTmp[color][q] =  fCell[color][q] - omega * (fCell[color][q] - fEq[color][q]);
            if (fTmp[color][q] < 0) fTmp[color][q] = 0;
        }
    } // end for 
    Cell result(fTmp);
    result.calcRho();
    return result;
}