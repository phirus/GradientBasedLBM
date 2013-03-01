#include<iostream>
#include<fstream>
#include"tests.h"

using namespace std;

//void Ausgabefile(Lattice);
//void Ausgabe(Lattice);

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();

}


//void Ausgabefile(Lattice D)
//{
//     ofstream rhofile;
//     int ysize = D.getYsize();
//     int xsize = D.getXsize();
//     Cell tmp;
//     double rho;
//
//     rhofile.open("rho.dat");
//     for(int j=0; j<ysize; j++) {
//          for(int i=0; i<xsize; i++) {
//               tmp = D.getCell(j,i);
//               tmp.calcRho(rho);
//               rhofile <<rho<<"\t";
//          }
//          rhofile<<endl;
//     }
//     rhofile.close();
//}
//
//void Ausgabe(Lattice D)
//{
//     int ysize = D.getYsize();
//     int xsize = D.getXsize();
//     Cell tmp;
//     double rho;
//
//     for(int j=0; j<ysize; j++) {
//          for(int i=0; i<xsize; i++) {
//               tmp = D.getCell(j,i);
//               tmp.calcRho(rho);
//               cout<<rho<<"\t";
//          }
//          cout<<endl;
//     }
//}
