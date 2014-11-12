#ifndef TESTS_H
#define TESTS_H

#include"gtest/gtest.h"
#include"../src/2D/Lattice2D.h"
#include"../src/2D/Vector2D.h"
#include"../src/BinaryIO.h"

#include"../src/3D/Definitions3D.h"
#include"../src/3D/Constants3D.h"

using namespace std;

TEST(BinaryIO,output){
    Lattice2D lattice(150,150);
    write_binary(lattice);
    Lattice2D vergleich;
    EXPECT_FALSE(read_binary(vergleich,"existiertnicht.txt"));
    EXPECT_TRUE(read_binary(vergleich));
    EXPECT_EQ(lattice,vergleich);
}

TEST(BinaryIO,paramLog){
    RelaxationPar rel = RelaxationPar(0.8,1.2,1.2);
    ParamSet params(1.1, 0.9, 1.1, 5, 2e-4, 2e-4, 1e-4, 1e-3,rel, 0.21, 0.11, 0.98);   
    EXPECT_NO_THROW(write_param_log(params));

    ParamSet read_in = read_paramset_file();
    EXPECT_EQ(params, read_in);
}

TEST(BinaryIO,queryTest){
    double value;
    EXPECT_FALSE( input_query("existiertnicht","test",value) );
    EXPECT_FALSE(input_query("queryTest","noflag",value));
    EXPECT_TRUE(input_query("queryTest","test",value));
    EXPECT_DOUBLE_EQ(13.4, value); 
}

TEST(BinaryIO,restart){
    Lattice2D lattice(150,150);

    Timetrack time(2e5,100,1000);
    time.timestep();

    Preprocess newProcess = read_preprocess_file("preprocessFile");
    write_restart_file(lattice, newProcess,time);
    
    Lattice2D vergleichL;
    Preprocess vergleichP;
    Timetrack vergleichT;
    EXPECT_TRUE(read_restart_file(vergleichL,vergleichP,vergleichT));
    EXPECT_EQ(lattice, vergleichL);
    EXPECT_EQ(newProcess, vergleichP);
    EXPECT_EQ(time, vergleichT);
}

TEST(Cell2D,constructor0)
{
    array2D f = {{1,0,0,0,0,0,0,0,0}};
    Cell2D cell;
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    EXPECT_EQ(f,cell.getF()[1]);
}

TEST(Cell2D,constructorRed1)
{
    array2D f = {{1,0,0,0,0,0,0,0,0}};
    Cell2D cell(1,0);
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell2D cell1(1,1);
    Cell2D cell2(1,100);
    EXPECT_EQ(f,cell1.getF()[0]);
    EXPECT_EQ(f,cell2.getF()[0]);
}

TEST(Cell2D,constructorBlue1)
{
    array2D f = {{1,0,0,0,0,0,0,0,0}};
    Cell2D cell(0,1);
    EXPECT_EQ(f,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell2D cell1(1,1);
    Cell2D cell2(100,1);
    EXPECT_EQ(f,cell1.getF()[1]);
    EXPECT_EQ(f,cell2.getF()[1]);
}

TEST(Cell2D,constructorRed2)
{
    array2D f = {{0,0,0,0,0,0,0,0,0}};
    Cell2D cell(0,0);
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell2D cell1(0,1);
    Cell2D cell2(0,100);
    EXPECT_EQ(f,cell1.getF()[0]);
    EXPECT_EQ(f,cell2.getF()[0]);
}

TEST(Cell2D,constructorBlue2)
{
    array2D f = {{0,0,0,0,0,0,0,0,0}};
    Cell2D cell(0,0);
    EXPECT_EQ(f,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell2D cell1(1,0);
    Cell2D cell2(100,0);
    EXPECT_EQ(f,cell1.getF()[1]);
    EXPECT_EQ(f,cell2.getF()[1]);
}

TEST(Cell2D,constructorSolid)
{
    Cell2D cell(0,0,true);
    EXPECT_EQ(true,cell.getIsSolid());
}

TEST(Cell2D,constructorIni)
{
    array2D f1 = {{4,2,5.4,0,1,12,6,7,8}};
    array2D f2 = {{1,2,3,4,5,6,7,8,9}};
    Cell2D cell(f1,f2);
    EXPECT_EQ(f1,cell.getF()[0]);
    EXPECT_EQ(f2,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
}

TEST(Cell2D,rho)
{
    array2D f1 = {{10,20,30,40,50,60,70,80,90}};
    array2D f2 = {{1,2,3,4,5,6,7,8,9}};
    Cell2D cell(f1,f2);
    EXPECT_EQ(0,cell.getRho()[0]);
    EXPECT_EQ(0,cell.getRho()[1]);
    cell.calcRho();
    EXPECT_EQ(450,cell.getRho()[0]);
    EXPECT_EQ(45,cell.getRho()[1]);
    cell.setIsSolid(true);
    cell.calcRho();
    EXPECT_EQ(0,cell.getRho()[0]);
    EXPECT_EQ(0,cell.getRho()[1]);
}

TEST(Cell2D,Psi1)
{
    array2D f1 = {{10,20,30,40,50,60,70,80,90}}; //450
    array2D f2 = {{1,2,3,4,10,6,7,8,9}}; //50
    Cell2D cell(f1,f2);
    EXPECT_DOUBLE_EQ(0,cell.calcPsi());
    cell.calcRho();
    EXPECT_DOUBLE_EQ(0.8, cell.calcPsi());
}

TEST(Cell2D,Psi2)
{
    array2D f1 = {{10,20,30,40,50,60,70,80,90}}; //450
    array2D f2 = {{0,0,0,0,0,0,0,0,}}; //0
    Cell2D cell(f1,f2);
    cell.calcRho();
    EXPECT_DOUBLE_EQ(1, cell.calcPsi());
    Cell2D cell2(f2,f1);
    cell2.calcRho();
    EXPECT_DOUBLE_EQ(-1, cell2.calcPsi());
}

TEST(Cell2D,Velo)
{
    array2D fr = {{2,3,0,5,0,0,0,0,0}};
    array2D fb = {{0,0,1,0,0,4,1,4,0}};

    double usqr;

    Cell2D cell(fr,fb);
    cell.calcRho();
    ColSet rho = cell.getRho();
    VeloSet2D v = cell.getU();

    Vector2D u = (v[0]*rho[0] + v[1]*rho[1]) * (1/sum(rho));

    usqr = u*u;

    EXPECT_DOUBLE_EQ(0.05,u.y);
    EXPECT_DOUBLE_EQ(-0.05,u.x);
    EXPECT_DOUBLE_EQ(0.005,(usqr));
}

TEST(Cell2D,equal){
    Cell2D one,two,three, four;
    EXPECT_EQ(one,two);
    three.setIsSolid(true);
    EXPECT_FALSE(one == three);
    array2D fr = {{2,3,0,5,0,0,0,0,0}};
    DistributionSetType2D f = four.getF();
    f[1] = fr;
    four.setF(f);
    EXPECT_FALSE(one == four);
}

TEST(Constants2D,BReis)
{
    EXPECT_DOUBLE_EQ(-4,B_2D[0]*27);
    EXPECT_DOUBLE_EQ(2,B_2D[1]*27);
    EXPECT_DOUBLE_EQ(5,B_2D[2]*108);
}

TEST(Constants2D,W)
{
    EXPECT_DOUBLE_EQ(4,WEIGHTS_2D.at(0)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_2D.at(1)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_2D.at(2)*36);
}

TEST(Constants3D,W)
{
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(0)*3);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(1)*18);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(2)*36);

    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(9)*18);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(14)*18);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(13)*36);
    EXPECT_DOUBLE_EQ(1,WEIGHTS_3D.at(12)*36);
}

TEST(Constants2D,Xi)
{
    EXPECT_DOUBLE_EQ(0,GRAD_WEIGHTS_2D.at(0));
    EXPECT_DOUBLE_EQ(32,GRAD_WEIGHTS_2D.at(1)*120);
    EXPECT_DOUBLE_EQ(12,GRAD_WEIGHTS_2D.at(2)*120);
    EXPECT_DOUBLE_EQ(1,GRAD_WEIGHTS_2D.at(9)*120);
}

TEST(Definitions2D,array2D_diff_add)
{
    const array2D one = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45}};
    const array2D two = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25}};
    const array2D vergleich = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20}};

    EXPECT_EQ (vergleich, array_diff_2D(one, two));
    EXPECT_EQ (one, array_add_2D(vergleich, two));
    EXPECT_EQ (one, array_add_2D(two, vergleich));
}

TEST(Definitions2D,distro_add_array2D)
{
    const array2D one = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20}};
    const array2D zeros = {{0,0,0,0,0,0,0,0,0}};
    DistributionSetType2D first = {{one, zeros}};

    const array2D plus = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25}};
    
    const array2D result = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45}};
    
    DistributionSetType2D vergleich = {{result, plus}};

    EXPECT_EQ(vergleich, distro_add_array_2D(first,plus));
}

TEST(Definitions2D,array2D_times)
{
    const array2D one = {{100, 5, 14, -2, 0, 0, 1, 916, 45}};
    const array2D vergleich = {{110, 5.5, 15.4, -2.2, 0, 0, 1.1, 1007.6, 49.5}};
    const array2D result = array_times_2D(one, 1.1);

    for(int i=0;i<9;i++){
        EXPECT_DOUBLE_EQ (vergleich[i], result[i]);
    }    
}

TEST(Definitions3D,array3D_diff_add)
{
    const array3D one = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45, 10, 11, 12, 13, 14, 15, 16, 17, 18,19}};
    const array3D two = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10}};
    const array3D vergleich = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20, -9, -7, -5, -3, -1, 1, 3, 5, 7, 9}};

    EXPECT_EQ (vergleich, array_diff_3D(one, two));
    EXPECT_EQ (one, array_add_3D(vergleich, two));
    EXPECT_EQ (one, array_add_3D(two, vergleich));
}

TEST(Definitions3D,distro_add_array3D)
{
    const array3D one = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}};
    const array3D zeros = {{0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0}};
    DistributionSetType3D first = {{one, zeros}};

    const array3D plus = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25,  19, 18, 17, 16, 15, 14, 13, 12, 11, 10}};
    
    const array3D result = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45, 29, 29, 29, 29, 29, 29, 29, 29, 29, 29}};
    
    DistributionSetType3D vergleich = {{result, plus}};

    EXPECT_EQ(vergleich, distro_add_array_3D(first,plus));
}

TEST(Definitions3D,array_times_3D)
{
    const array3D one = {{100, 5, 14, -2, 0, 0, 1, 916, 45, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19}};
    const array3D vergleich = {{110, 5.5, 15.4, -2.2, 0, 0, 1.1, 1007.6, 49.5, 11, 12.1, 13.2, 14.3, 15.4, 16.5, 17.6, 18.7, 19.8, 20.9}};
    const array3D result = array_times_3D(one, 1.1);

    for(int i=0;i<19;i++){
        EXPECT_DOUBLE_EQ (vergleich[i], result[i]);
    }    
}

TEST(Lattice2D,constructor)
{
    const Lattice2D lattice;
    EXPECT_EQ(10, lattice.getSize()[1]);
    EXPECT_EQ(10, lattice.getSize()[0]);
    EXPECT_EQ(1,lattice.getF(1,1)[0][0]);
}

TEST(Lattice2D,stream)
{
    Lattice2D lattice(3,3,0,0);

    array2D f1 = {{0,1,0,0,0,0,0,0,0}};
    array2D f2 = {{0,0,1,0,0,0,0,0,0}};
    array2D f3 = {{0,0,0,1,0,0,0,0,0}};
    array2D f4 = {{0,0,0,0,1,0,0,0,0}};
    array2D f5 = {{0,0,0,0,0,1,0,0,0}};
    array2D f6 = {{0,0,0,0,0,0,1,0,0}};
    array2D f7 = {{0,0,0,0,0,0,0,1,0}};
    array2D f8 = {{0,0,0,0,0,0,0,0,1}};

    lattice.setF(0,0,0,f2);
    lattice.setF(1,0,0,f3);
    lattice.setF(2,0,0,f4);
    lattice.setF(0,1,0,f1);
    lattice.setF(2,1,0,f5);
    lattice.setF(0,2,0,f8);
    lattice.setF(1,2,0,f7);
    lattice.setF(2,2,0,f6);

    array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(1,1,1,fcenter);

    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(1,1)[0]);

    EXPECT_EQ(f6,lattice.getF(0,0)[1]);
    EXPECT_EQ(f7,lattice.getF(1,0)[1]);
    EXPECT_EQ(f8,lattice.getF(2,0)[1]);
    EXPECT_EQ(f5,lattice.getF(0,1)[1]);
    EXPECT_EQ(f1,lattice.getF(2,1)[1]);
    EXPECT_EQ(f4,lattice.getF(0,2)[1]);
    EXPECT_EQ(f3,lattice.getF(1,2)[1]);
    EXPECT_EQ(f2,lattice.getF(2,2)[1]);
}

TEST(Lattice2D,bounceSmall)
{
    Lattice2D lattice(3,3,0,0);
    lattice.closedBox();

    array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    array2D fzero = {{0,0,0,0,0,0,0,0,0}};

    lattice.setF(1,1,0,fcenter);
    lattice.setF(1,1,1,fcenter);

    lattice.streamAll();

    EXPECT_EQ(fzero,lattice.getF(0,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(0,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(0,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(1,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(1,2)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,2)[1]);

    EXPECT_EQ(fzero,lattice.getF(0,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(0,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(0,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(1,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(1,2)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,0)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,1)[1]);
    EXPECT_EQ(fzero,lattice.getF(2,2)[1]);

    EXPECT_EQ(fcenter,lattice.getF(1,1)[0]);
    EXPECT_EQ(fcenter,lattice.getF(1,1)[1]);
}

TEST(Lattice2D,bounceClosed)
{
    Lattice2D lattice(5,5,0,0);
    lattice.closedBox();
    lattice.setCell(2,2,Cell2D(0,0,true));
    Cell2D tmp = lattice.getCell(2,2);
    EXPECT_EQ(true,tmp.getIsSolid());

    const array2D f1 = {{0,1,0,0,0,0,0,0,0}};
    const array2D f2 = {{0,0,1,0,0,0,0,0,0}};
    const array2D f3 = {{0,0,0,1,0,0,0,0,0}};
    const array2D f4 = {{0,0,0,0,1,0,0,0,0}};
    const array2D f5 = {{0,0,0,0,0,1,0,0,0}};
    const array2D f6 = {{0,0,0,0,0,0,1,0,0}};
    const array2D f7 = {{0,0,0,0,0,0,0,1,0}};
    const array2D f8 = {{0,0,0,0,0,0,0,0,1}};

    lattice.setF(1,1,0,f2);
    lattice.setF(2,1,0,f3);
    lattice.setF(3,1,0,f4);
    lattice.setF(1,2,0,f1);
    lattice.setF(3,2,0,f5);
    lattice.setF(1,3,0,f8);
    lattice.setF(2,3,0,f7);
    lattice.setF(3,3,0,f6);

    lattice.setF(1,1,1,f2);
    lattice.setF(2,1,1,f3);
    lattice.setF(3,1,1,f4);
    lattice.setF(1,2,1,f1);
    lattice.setF(3,2,1,f5);
    lattice.setF(1,3,1,f8);
    lattice.setF(2,3,1,f7);
    lattice.setF(3,3,1,f6);

    lattice.streamAll();

    EXPECT_EQ(f6,lattice.getF(1,1)[0]);
    EXPECT_EQ(f7,lattice.getF(2,1)[0]);
    EXPECT_EQ(f8,lattice.getF(3,1)[0]);
    EXPECT_EQ(f5,lattice.getF(1,2)[0]);
    EXPECT_EQ(f1,lattice.getF(3,2)[0]);
    EXPECT_EQ(f4,lattice.getF(1,3)[0]);
    EXPECT_EQ(f3,lattice.getF(2,3)[0]);
    EXPECT_EQ(f2,lattice.getF(3,3)[0]);

    EXPECT_EQ(f6,lattice.getF(1,1)[1]);
    EXPECT_EQ(f7,lattice.getF(2,1)[1]);
    EXPECT_EQ(f8,lattice.getF(3,1)[1]);
    EXPECT_EQ(f5,lattice.getF(1,2)[1]);
    EXPECT_EQ(f1,lattice.getF(3,2)[1]);
    EXPECT_EQ(f4,lattice.getF(1,3)[1]);
    EXPECT_EQ(f3,lattice.getF(2,3)[1]);
    EXPECT_EQ(f2,lattice.getF(3,3)[1]);
}

TEST(Lattice2D,bounceClosed2)
{
    Lattice2D lattice(5,5,0,0);
    lattice.closedBox();
    array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,0,fcenter);
    lattice.setF(2,2,1,fcenter);
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(2,2)[0]);
    EXPECT_EQ(fcenter,lattice.getF(2,2)[1]);
}

TEST(Lattice2D,periodic)
{
    Lattice2D lattice(5,5,0,0);
    array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,0,fcenter);
    lattice.setF(2,2,1,fcenter);

    // 5 mal
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(2,2)[0]);
    EXPECT_EQ(fcenter,lattice.getF(2,2)[1]);
}

TEST(Lattice2D,streamRho)
{
    Cell2D tmp;

    Lattice2D lattice(5,5,0,0);
    lattice.closedBox();

    array2D fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,0,fcenter);
    lattice.setF(2,2,1,fcenter);
    lattice.streamAll();

    tmp = lattice.getCell(2,2);
    EXPECT_EQ(0,tmp.getRho()[1]);
    tmp = lattice.getCell(1,2);
    EXPECT_EQ(1,tmp.getRho()[1]);
    tmp = lattice.getCell(3,2);
    EXPECT_EQ(1,tmp.getRho()[1]);

    lattice.streamAll();
    tmp = lattice.getCell(1,2);
    EXPECT_EQ(1,tmp.getRho()[1]);
    tmp = lattice.getCell(3,2);
    EXPECT_EQ(1,tmp.getRho()[1]);
}

TEST(Lattice2D,collideSingle)
{
    Lattice2D lattice(5,5,1,1);
    lattice.closedBox();

    lattice.collideAll();
    Cell2D cell = lattice.getCell(2,2);
    double rho,usqr;
    cell.calcRho();
    rho = sum(cell.getRho());
    VeloSet2D u = cell.getU();
    usqr = u[0]*u[0];

    EXPECT_DOUBLE_EQ(2.0,rho);
    EXPECT_NEAR(0,usqr,1e-10);
}

TEST(Lattice2D, directions)
{
    Lattice2D lattice(5,5);
    direction2D dir = lattice.directions(0,0);

    boost::array<int,13> x = {{0,1,1,0,4,4,4,0,1,2,0,3,0}};
    boost::array<int,13> y = {{0,0,1,1,1,0,4,4,4,0,2,0,3}};
    for(int q=0; q<13;q++){
    EXPECT_EQ(y[q],dir[q].y);
    EXPECT_EQ(x[q],dir[q].x);
    }
}

TEST(Lattice2D, Gradient)
{
    Lattice2D lattice(5,5,0,0);
    array2D f = {{1,0,0,0,0,0,0,0,0}};

    for (int i=0; i<5; i++)
    {
        lattice.setF(0,i,1,f);
        lattice.setF(1,i,1,f);
        lattice.setF(2,i,1,f);
        lattice.setF(3,i,0,f);
        lattice.setF(4,i,0,f);
    }
    Cell2D cell;
    for (int i = 0; i<5; i++)
    {
        for (int j = 0; j<5; j++)
        {
            cell = lattice.getCell(j,i);
            cell.calcRho();
            lattice.setCell(j,i,cell);
        }

    }
    Vector2D grad = lattice.getGradient(2,2);
    double abs = grad.Abs();
    grad.x /= abs;
    grad.y /= abs;

    EXPECT_DOUBLE_EQ(1,grad.x); // y-Richtung
    EXPECT_DOUBLE_EQ(0,grad.y); // x-Richung
}

TEST(Lattice2D,collisionBalanceAll)
{
    Lattice2D lattice(5,5,1,1);
    lattice.closedBox();
    double mass, momentum;

    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(18,mass);
    EXPECT_DOUBLE_EQ(0,momentum);

    lattice.collideAll();
    lattice.streamAll();
    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(18,mass);
    EXPECT_NEAR(0,momentum,1e-10);

    lattice.collideAll();
    lattice.streamAll();

    lattice.balance(mass, momentum);
    EXPECT_DOUBLE_EQ(18,mass);
    EXPECT_NEAR(0,momentum,1e-10);
}

TEST(Lattice2D, copy_constr){
    Lattice2D lBig(100,20);
    Lattice2D tmp(lBig);
    EXPECT_EQ(lBig,tmp);
}

TEST(Lattice2D, assign){ 
    Lattice2D lBig(100,20);
    Lattice2D tmp;
    tmp = lBig;
    EXPECT_EQ(lBig,tmp);
}

TEST(Matrix2D,trafo){
    const array2D verteilung = {{1,2,3,4,5,6,7,8,9}};
    const array2D vergleich = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array2D trafo = TRAFO_MATRIX * verteilung;

    EXPECT_EQ(vergleich, trafo);
}

TEST(Matrix2D,backtrafo){
    const array2D vergleich = {{1,2,3,4,5,6,7,8,9}};
    const array2D verteilung = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array2D backtrafo = INV_TRAFO_MATRIX * verteilung;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],backtrafo[i]);
    }
}

TEST(Matrix2D,multiply){
    const Matrix2D S(RelaxationPar(1,10,100));
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    // const array2D vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};
    const array2D vergleich = {{ 1, 2, 30, 4, 500, 6, 700, 8, 9}};

    array2D test = S*f;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],test[i])<<"i = "<<i ;
    }
}

TEST(Matrix2D,multiply_linewise){
    const Matrix2D S(RelaxationPar(1,10,100));
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    // const array2D vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};
    const array2D vergleich = {{ 1, 2, 30, 4, 500, 6, 700, 8, 9}};

    array2D test;

    for(int i = 0; i<9;i++)
    {
        test[i] = S.linewise(f,i);
        EXPECT_DOUBLE_EQ(vergleich[i],test[i])<<"i = "<<i ;
    }
}

TEST(Matrix2D,identity){
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    const array2D f0 = {{0,0,0,0,0,0,0,0,0}};

    const Matrix2D Identity = Matrix2D(true);
    const Matrix2D Zeros = Matrix2D();
    array2D test;

    EXPECT_EQ(f, Identity * f);
    EXPECT_EQ(f0, Zeros * f);
    for(int i = 0; i<9;i++)
    {
        test[i] = Identity.linewise(f,i);
    }
    EXPECT_EQ(f, test);
}

TEST(Matrix2D,plus_times){
    const Matrix2D Identity = Matrix2D(true);
    
    EXPECT_EQ(Identity+Identity, Identity*2);

    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX+TRAFO_MATRIX, TRAFO_MATRIX*3);
    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX, (TRAFO_MATRIX*3)-TRAFO_MATRIX);
}

TEST(MRT,trafo){
    /// testet ob die Differenz im Geschw.-Raum gleich der Rücktransformierten Differenz im moment-Raum ist
    ParamSet param;
    DistributionSetType2D phi = param.getPhi();
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    Cell2D testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();

    VeloSet2D u = testCell.getU();
    DistributionSetType2D fEq  = eqDistro(rho,u,phi);
    DistributionSetType2D vergleich;
    vergleich[0] = array_diff_2D(f, fEq[0]);
    vergleich[1] = array_diff_2D(f, fEq[1]);
    array2D m = TRAFO_MATRIX * f;
    DistributionSetType2D mEq;
    mEq[0] = TRAFO_MATRIX * fEq[0];
    mEq[1] = TRAFO_MATRIX * fEq[1];
    DistributionSetType2D transformed;
    transformed[0] = INV_TRAFO_MATRIX * array_diff_2D(m,mEq[0]);
    transformed[1] = INV_TRAFO_MATRIX * array_diff_2D(m,mEq[1]);

    for(int i = 0; i<9;i++)
    {
        EXPECT_NEAR(vergleich[0][i],transformed[0][i],1e-10) ;
        EXPECT_NEAR(vergleich[1][i],transformed[1][i],1e-10) ;
    }
}

TEST(MRT,mass){
    ParamSet param;
    DistributionSetType2D phi = param.getPhi();
    const array2D f = {{1,2,3,4,5,6,7,8,9}};
    Cell2D testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();
    VeloSet2D u = testCell.getU();
    DistributionSetType2D fEq  = eqDistro(rho,u,phi);
    array2D m = TRAFO_MATRIX * f;
    array2D mEq = TRAFO_MATRIX * fEq[0];

    EXPECT_DOUBLE_EQ(m[0],mEq[0]);
    EXPECT_DOUBLE_EQ(m[3],mEq[3]);
    EXPECT_DOUBLE_EQ(m[5],mEq[5]);
}

TEST(ParamSet,Phi)
{
    ParamSet param(1,1,1,1000);
    DistributionSetType2D phi;
    phi = param.getPhi();

    array2D phiR = {{0.9992, 1.6e-4, 4e-5, 1.6e-4, 4e-5, 1.6e-4, 4e-5, 1.6e-4, 4e-5}};
    array2D phiB = {{0.2,0.16,0.04,0.16,0.04,0.16,0.04,0.16,0.04}};

    for(int q=0;q<9;q++)
    {
        EXPECT_NEAR(phiR[q], phi.at(0)[q],1e-10);
    }
    EXPECT_EQ(phiB, phi.at(1));
}

TEST(ParamSet,inter)
{
    ParamSet param;
    param.setOmega(1.1,0.9,0.1);
    Interpol inter = param.getInter();

    EXPECT_DOUBLE_EQ(0.99,inter.chi);
    EXPECT_DOUBLE_EQ(2.2,inter.eta);
    EXPECT_DOUBLE_EQ(-11,inter.kappa);
    EXPECT_NEAR(1.8,inter.lambda,1e-10);
    EXPECT_NEAR(9,inter.ny,1e-10);
}

TEST(ParamSet, Ak)
{
    double omega(1.2), sigma(1.5e-4);
    ParamSet param(1, 1, 1, 2, sigma);
    ColSet A_k = param.getAk(omega);
    EXPECT_DOUBLE_EQ(sigma, param.getSigma());
    EXPECT_DOUBLE_EQ(0.0004049999999999999, A_k[0]);
    EXPECT_DOUBLE_EQ(0.0004049999999999999, A_k[1]);
}

TEST(ParamSet, equal){
    ParamSet one, two, three, four;
    EXPECT_TRUE(one == two);
    three.setBeta(0);
    EXPECT_FALSE(one == three);
    four.setRelaxation(1,1,1);
    EXPECT_EQ(one,four);
}

TEST(Preprocess,constr){
    Preprocess newProcess;

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(10,newProcess.getReynoldsMax());
    EXPECT_DOUBLE_EQ(100,newProcess.getMorton());
    EXPECT_DOUBLE_EQ(10,newProcess.getEotvos());
    EXPECT_EQ(30,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1,newProcess.getRhoL());
    EXPECT_DOUBLE_EQ(2,newProcess.getGamma());
    EXPECT_DOUBLE_EQ(2,newProcess.getMuRatio());

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(1,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(1,newProcess.getTimestep());

    EXPECT_DOUBLE_EQ(1.0196152422706632,newProcess.getTau());
    EXPECT_DOUBLE_EQ(0.5,newProcess.getDelRho());
    EXPECT_DOUBLE_EQ(0.5773502691896258,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(0.17320508075688779,newProcess.getNu());
    EXPECT_DOUBLE_EQ(0.000316227766016838,newProcess.getSigma());
    EXPECT_DOUBLE_EQ(7.0272836892630665e-6,newProcess.getG());

    EXPECT_EQ(120,newProcess.getWidth());
    EXPECT_EQ(360,newProcess.getHeight());
}

TEST(Preprocess,FileInput){
    Preprocess newProcess = read_preprocess_file("preprocessFile");

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(15,newProcess.getReynoldsMax());
    EXPECT_DOUBLE_EQ(10,newProcess.getMorton());
    EXPECT_DOUBLE_EQ(1,newProcess.getEotvos());
    EXPECT_EQ(35,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1,newProcess.getRhoL());
    EXPECT_DOUBLE_EQ(3,newProcess.getGamma());
    EXPECT_DOUBLE_EQ(1.8,newProcess.getMuRatio());

    // test the stored parameters
    EXPECT_DOUBLE_EQ(1,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(1,newProcess.getTimestep());
    EXPECT_EQ(150,newProcess.getWidth());
    EXPECT_EQ(400,newProcess.getHeight());

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(0.9041451884327381,newProcess.getTau());
    EXPECT_DOUBLE_EQ(0.6666666666666667,newProcess.getDelRho());
    EXPECT_DOUBLE_EQ(0.5773502691896258,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(0.13471506281091272,newProcess.getNu());
    EXPECT_DOUBLE_EQ(0.00016396995274947155,newProcess.getSigma());
    EXPECT_DOUBLE_EQ(2.0077953397894473e-7,newProcess.getG());
}

TEST(Preprocess,ParameterCheck){
    Preprocess newProcess = read_preprocess_file("realParameters");
    ParamSet params = newProcess.getParamSet();

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(75,newProcess.getReynoldsMax());
    EXPECT_NEAR(0.01,newProcess.getMorton(),1e-10);
    EXPECT_DOUBLE_EQ(10,newProcess.getEotvos());

    EXPECT_DOUBLE_EQ(1.6558400984359005, params.getOmegaRed());
    EXPECT_DOUBLE_EQ(1.6558400984359005, params.getOmegaBlue());
    EXPECT_DOUBLE_EQ(1,params.getRhoR());
    EXPECT_DOUBLE_EQ(1.0295861249616975,params.getGamma());
    EXPECT_DOUBLE_EQ(0.2,params.getAlpha());
    EXPECT_DOUBLE_EQ(0.1,params.getInterfaceThickness());
    EXPECT_DOUBLE_EQ(0.99,params.getBeta());
    // EXPECT_DOUBLE_EQ(0.00084327404271156506,params.getSigma());
    // EXPECT_DOUBLE_EQ(0.00014491665424627533,params.getG());
     }

TEST(timetrack,FileInput){
    Preprocess newProcess = read_preprocess_file("preprocessFile");
    Timetrack newTimetrack = read_timetrack_file("preprocessFile");
     // test the given parameters 
    EXPECT_DOUBLE_EQ(1e6,newTimetrack.getMaxCount());
    EXPECT_DOUBLE_EQ(1000,newTimetrack.getOutputInt());
    EXPECT_DOUBLE_EQ(10000,newTimetrack.getRestartInt());
}

TEST(Vector2D,scalar){
    Vector2D v0, v1(1,2), v2(3,4);
    EXPECT_DOUBLE_EQ(0, v0*v1);
    EXPECT_DOUBLE_EQ(0, v2*v0);
    EXPECT_DOUBLE_EQ(11, v1*v2);
    EXPECT_DOUBLE_EQ(11, v2*v1);
}

TEST(Vector2D,angle){
    Vector2D g(1,1);

    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[0]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION_2D[1]));
    EXPECT_DOUBLE_EQ(1, g.Angle(DIRECTION_2D[2]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION_2D[3]));
    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[4]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION_2D[5]));
    EXPECT_DOUBLE_EQ(-1, g.Angle(DIRECTION_2D[6]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION_2D[7]));
    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[8]));

    Vector2D g1(1e-10, -1e-10), g2(1e-6, -1e-6);
    EXPECT_DOUBLE_EQ(1,g1.Angle(g2));
    EXPECT_DOUBLE_EQ(0,g1.Angle(DIRECTION_2D[0]));
}

TEST(Vector3D,scalar){
    Vector3D v0, v1(1,2,3), v2(4,5,6), v3(-1,-20,100);
    EXPECT_DOUBLE_EQ(0, v0*v1);
    EXPECT_DOUBLE_EQ(0, v2*v0);
    EXPECT_DOUBLE_EQ(32, v1*v2);
    EXPECT_DOUBLE_EQ(32, v2*v1);
    EXPECT_DOUBLE_EQ(259, v1*v3);
}

// TEST(Vector2D,angle){
//     Vector2D g(1,1);

//     EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[0]));
//     EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION_2D[1]));
//     EXPECT_DOUBLE_EQ(1, g.Angle(DIRECTION_2D[2]));
//     EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION_2D[3]));
//     EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[4]));
//     EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION_2D[5]));
//     EXPECT_DOUBLE_EQ(-1, g.Angle(DIRECTION_2D[6]));
//     EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION_2D[7]));
//     EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION_2D[8]));

//     Vector2D g1(1e-10, -1e-10), g2(1e-6, -1e-6);
//     EXPECT_DOUBLE_EQ(1,g1.Angle(g2));
//     EXPECT_DOUBLE_EQ(0,g1.Angle(DIRECTION_2D[0]));
// }

#endif