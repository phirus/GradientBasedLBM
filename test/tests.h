#ifndef TESTS_H
#define TESTS_H

#include"gtest/gtest.h"
#include"../src/Lattice.h"
#include"../src/Vector.h"
#include"../src/BinaryIO.h"

using namespace std;

TEST(BinaryIO,output){
    Lattice lattice(150,150);
    write_binary(lattice);
    Lattice vergleich;
    EXPECT_FALSE(read_binary(vergleich,"existiertnicht.txt"));
    EXPECT_TRUE(read_binary(vergleich));
    EXPECT_EQ(lattice,vergleich);
}

TEST(BinaryIO,paramLog){
    RelaxationPar rel = RelaxationPar(0.8,1.2,1.2);
    ParamSet params(1.1, 0.9, 1.1, 5, 2e-4, 2e-4, 10, 1e-4, rel, 0.21, 0.11, 0.98);   
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
    Lattice lattice(150,150);

    Timetrack time(1e-2, 1.05, 2e5,100,1000);
    time.timestep();
    time.timestep();
    time.timestep();
    time.refine();
    time.timestep();
    time.timestep();

    Preprocess newProcess = read_preprocess_file("preprocessFile");
    write_restart_file(lattice, newProcess,time);
    
    Lattice vergleichL;
    Preprocess vergleichP;
    Timetrack vergleichT;
    EXPECT_TRUE(read_restart_file(vergleichL,vergleichP,vergleichT));
    EXPECT_EQ(lattice, vergleichL);
    EXPECT_EQ(newProcess, vergleichP);
    EXPECT_EQ(time, vergleichT);
}

TEST(Cell,constructor0)
{
    array f = {{1,0,0,0,0,0,0,0,0}};
    Cell cell;
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    EXPECT_EQ(f,cell.getF()[1]);
}

TEST(Cell,constructorRed1)
{
    array f = {{1,0,0,0,0,0,0,0,0}};
    Cell cell(1,0);
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell cell1(1,1);
    Cell cell2(1,100);
    EXPECT_EQ(f,cell1.getF()[0]);
    EXPECT_EQ(f,cell2.getF()[0]);
}

TEST(Cell,constructorBlue1)
{
    array f = {{1,0,0,0,0,0,0,0,0}};
    Cell cell(0,1);
    EXPECT_EQ(f,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell cell1(1,1);
    Cell cell2(100,1);
    EXPECT_EQ(f,cell1.getF()[1]);
    EXPECT_EQ(f,cell2.getF()[1]);
}

TEST(Cell,constructorRed2)
{
    array f = {{0,0,0,0,0,0,0,0,0}};
    Cell cell(0,0);
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell cell1(0,1);
    Cell cell2(0,100);
    EXPECT_EQ(f,cell1.getF()[0]);
    EXPECT_EQ(f,cell2.getF()[0]);
}

TEST(Cell,constructorBlue2)
{
    array f = {{0,0,0,0,0,0,0,0,0}};
    Cell cell(0,0);
    EXPECT_EQ(f,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
    Cell cell1(1,0);
    Cell cell2(100,0);
    EXPECT_EQ(f,cell1.getF()[1]);
    EXPECT_EQ(f,cell2.getF()[1]);
}

TEST(Cell,constructorSolid)
{
    Cell cell(0,0,true);
    EXPECT_EQ(true,cell.getIsSolid());
}

TEST(Cell,constructorIni)
{
    array f1 = {{4,2,5.4,0,1,12,6,7,8}};
    array f2 = {{1,2,3,4,5,6,7,8,9}};
    Cell cell(f1,f2);
    EXPECT_EQ(f1,cell.getF()[0]);
    EXPECT_EQ(f2,cell.getF()[1]);
    EXPECT_FALSE(cell.getIsSolid());
}

TEST(Cell,rho)
{
    array f1 = {{10,20,30,40,50,60,70,80,90}};
    array f2 = {{1,2,3,4,5,6,7,8,9}};
    Cell cell(f1,f2);
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

TEST(Cell,Psi1)
{
    array f1 = {{10,20,30,40,50,60,70,80,90}}; //450
    array f2 = {{1,2,3,4,10,6,7,8,9}}; //50
    Cell cell(f1,f2);
    EXPECT_DOUBLE_EQ(0,cell.calcPsi());
    cell.calcRho();
    EXPECT_DOUBLE_EQ(0.8, cell.calcPsi());
}

TEST(Cell,Psi2)
{
    array f1 = {{10,20,30,40,50,60,70,80,90}}; //450
    array f2 = {{0,0,0,0,0,0,0,0,}}; //0
    Cell cell(f1,f2);
    cell.calcRho();
    EXPECT_DOUBLE_EQ(1, cell.calcPsi());
    Cell cell2(f2,f1);
    cell2.calcRho();
    EXPECT_DOUBLE_EQ(-1, cell2.calcPsi());
}

TEST(Cell,Velo)
{
    array fr = {{2,3,0,5,0,0,0,0,0}};
    array fb = {{0,0,1,0,0,4,1,4,0}};

    double usqr;

    Cell cell(fr,fb);
    cell.calcRho();
    ColSet rho = cell.getRho();
    VeloSet v = cell.getU();

    Vector u = (v[0]*rho[0] + v[1]*rho[1]) * (1/sum(rho));

    usqr = u*u;

    EXPECT_DOUBLE_EQ(0.05,u.y);
    EXPECT_DOUBLE_EQ(-0.05,u.x);
    EXPECT_DOUBLE_EQ(0.005,(usqr));
}

TEST(Cell,equal){
    Cell one,two,three, four;
    EXPECT_EQ(one,two);
    three.setIsSolid(true);
    EXPECT_FALSE(one == three);
    array fr = {{2,3,0,5,0,0,0,0,0}};
    DistributionSetType f = four.getF();
    f[1] = fr;
    four.setF(f);
    EXPECT_FALSE(one == four);
}

TEST(Constants,BReis)
{
    EXPECT_DOUBLE_EQ(-4,B[0]*27);
    EXPECT_DOUBLE_EQ(2,B[1]*27);
    EXPECT_DOUBLE_EQ(5,B[2]*108);
}

TEST(Constants,W)
{
    EXPECT_DOUBLE_EQ(4,WEIGHTS.at(0)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS.at(1)*9);
    EXPECT_DOUBLE_EQ(1,WEIGHTS.at(2)*36);
}

TEST(Constants,Xi)
{
    EXPECT_DOUBLE_EQ(0,GRAD_WEIGHTS.at(0));
    EXPECT_DOUBLE_EQ(32,GRAD_WEIGHTS.at(1)*120);
    EXPECT_DOUBLE_EQ(12,GRAD_WEIGHTS.at(2)*120);
    EXPECT_DOUBLE_EQ(1,GRAD_WEIGHTS.at(9)*120);
}

TEST(Definitions,array_diff_add)
{
    const array one = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45}};
    const array two = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25}};
    const array vergleich = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20}};

    EXPECT_EQ (vergleich, array_diff(one, two));
    EXPECT_EQ (one, array_add(vergleich, two));
    EXPECT_EQ (one, array_add(two, vergleich));
}

TEST(Definitions,distro_add_array)
{
    const array one = {{95, 4.9, -0.5, 0, -1.8, -100, 0.5, 0, 20}};
    const array zeros = {{0,0,0,0,0,0,0,0,0}};
    DistributionSetType first = {{one, zeros}};

    const array plus = {{5,   0.6, 1.9, -2, 1.8, 100, 0.5, 91.6, 25}};
    
    const array result = {{100, 5.5, 1.4, -2, 0, 0, 1, 91.6, 45}};
    
    DistributionSetType vergleich = {{result, plus}};

    EXPECT_EQ(vergleich, distro_add_array(first,plus));
}

TEST(Definitions,array_times)
{
    const array one = {{100, 5, 14, -2, 0, 0, 1, 916, 45}};
    const array vergleich = {{110, 5.5, 15.4, -2.2, 0, 0, 1.1, 1007.6, 49.5}};
    const array result = array_times(one, 1.1);

    for(int i=0;i<9;i++){
        EXPECT_DOUBLE_EQ (vergleich[i], result[i]);
    }    
}

TEST(Lattice,constructor)
{
    const Lattice lattice;
    EXPECT_EQ(10, lattice.getSize()[1]);
    EXPECT_EQ(10, lattice.getSize()[0]);
    EXPECT_EQ(1,lattice.getF(1,1)[0][0]);
}

TEST(Lattice,stream)
{
    Lattice lattice(3,3,0,0);

    array f1 = {{0,1,0,0,0,0,0,0,0}};
    array f2 = {{0,0,1,0,0,0,0,0,0}};
    array f3 = {{0,0,0,1,0,0,0,0,0}};
    array f4 = {{0,0,0,0,1,0,0,0,0}};
    array f5 = {{0,0,0,0,0,1,0,0,0}};
    array f6 = {{0,0,0,0,0,0,1,0,0}};
    array f7 = {{0,0,0,0,0,0,0,1,0}};
    array f8 = {{0,0,0,0,0,0,0,0,1}};

    lattice.setF(0,0,0,f2);
    lattice.setF(1,0,0,f3);
    lattice.setF(2,0,0,f4);
    lattice.setF(0,1,0,f1);
    lattice.setF(2,1,0,f5);
    lattice.setF(0,2,0,f8);
    lattice.setF(1,2,0,f7);
    lattice.setF(2,2,0,f6);

    array fcenter = {{0,1,1,1,1,1,1,1,1}};
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

TEST(Lattice,bounceSmall)
{
    Lattice lattice(3,3,0,0);
    lattice.closedBox();

    array fcenter = {{0,1,1,1,1,1,1,1,1}};
    array fzero = {{0,0,0,0,0,0,0,0,0}};

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

TEST(Lattice,bounceClosed)
{
    Lattice lattice(5,5,0,0);
    lattice.closedBox();
    lattice.setCell(2,2,Cell(0,0,true));
    Cell tmp = lattice.getCell(2,2);
    EXPECT_EQ(true,tmp.getIsSolid());

    const array f1 = {{0,1,0,0,0,0,0,0,0}};
    const array f2 = {{0,0,1,0,0,0,0,0,0}};
    const array f3 = {{0,0,0,1,0,0,0,0,0}};
    const array f4 = {{0,0,0,0,1,0,0,0,0}};
    const array f5 = {{0,0,0,0,0,1,0,0,0}};
    const array f6 = {{0,0,0,0,0,0,1,0,0}};
    const array f7 = {{0,0,0,0,0,0,0,1,0}};
    const array f8 = {{0,0,0,0,0,0,0,0,1}};

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

TEST(Lattice,bounceClosed2)
{
    Lattice lattice(5,5,0,0);
    lattice.closedBox();
    array fcenter = {{0,1,1,1,1,1,1,1,1}};
    lattice.setF(2,2,0,fcenter);
    lattice.setF(2,2,1,fcenter);
    lattice.streamAll();
    lattice.streamAll();
    lattice.streamAll();

    EXPECT_EQ(fcenter,lattice.getF(2,2)[0]);
    EXPECT_EQ(fcenter,lattice.getF(2,2)[1]);
}

TEST(Lattice,periodic)
{
    Lattice lattice(5,5,0,0);
    array fcenter = {{0,1,1,1,1,1,1,1,1}};
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

TEST(Lattice,streamRho)
{
    Cell tmp;

    Lattice lattice(5,5,0,0);
    lattice.closedBox();

    array fcenter = {{0,1,1,1,1,1,1,1,1}};
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

TEST(Lattice,collideSingle)
{
    Lattice lattice(5,5,1,1);
    lattice.closedBox();

    lattice.collideAll();
    Cell cell = lattice.getCell(2,2);
    double rho,usqr;
    cell.calcRho();
    rho = sum(cell.getRho());
    VeloSet u = cell.getU();
    usqr = u[0]*u[0];

    EXPECT_DOUBLE_EQ(2.0,rho);
    EXPECT_NEAR(0,usqr,1e-10);
}

TEST(Lattice, directions)
{
    Lattice lattice(5,5);
    direction dir = lattice.directions(0,0);

    boost::array<int,13> x = {{0,1,1,0,4,4,4,0,1,2,0,3,0}};
    boost::array<int,13> y = {{0,0,1,1,1,0,4,4,4,0,2,0,3}};
    for(int q=0; q<13;q++){
    EXPECT_EQ(y[q],dir[q].y);
    EXPECT_EQ(x[q],dir[q].x);
    }
}

TEST(Lattice, Gradient)
{
    Lattice lattice(5,5,0,0);
    array f = {{1,0,0,0,0,0,0,0,0}};

    for (int i=0; i<5; i++)
    {
        lattice.setF(0,i,1,f);
        lattice.setF(1,i,1,f);
        lattice.setF(2,i,1,f);
        lattice.setF(3,i,0,f);
        lattice.setF(4,i,0,f);
    }
    Cell cell;
    for (int i = 0; i<5; i++)
    {
        for (int j = 0; j<5; j++)
        {
            cell = lattice.getCell(j,i);
            cell.calcRho();
            lattice.setCell(j,i,cell);
        }

    }
    Vector grad = lattice.getGradient(2,2);
    double abs = grad.Abs();
    grad.x /= abs;
    grad.y /= abs;

    EXPECT_DOUBLE_EQ(1,grad.x); // y-Richtung
    EXPECT_DOUBLE_EQ(0,grad.y); // x-Richung
}

TEST(Lattice,collisionBalanceAll)
{
    Lattice lattice(5,5,1,1);
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

TEST(Lattice, copy_constr){
    Lattice lBig(100,20);
    Lattice tmp(lBig);
    EXPECT_EQ(lBig,tmp);
}

TEST(Lattice, assign){ 
    Lattice lBig(100,20);
    Lattice tmp;
    tmp = lBig;
    EXPECT_EQ(lBig,tmp);
}

TEST(Matrix,trafo){
    const array verteilung = {{1,2,3,4,5,6,7,8,9}};
    const array vergleich = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array trafo = TRAFO_MATRIX * verteilung;

    EXPECT_EQ(vergleich, trafo);
}

TEST(Matrix,backtrafo){
    const array vergleich = {{1,2,3,4,5,6,7,8,9}};
    const array verteilung = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array backtrafo = INV_TRAFO_MATRIX * verteilung;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],backtrafo[i]);
    }
}

TEST(Matrix,multiply){
    const Matrix S(RelaxationPar(1,10,100));
    const array f = {{1,2,3,4,5,6,7,8,9}};
    // const array vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};
    const array vergleich = {{ 1, 2, 30, 4, 500, 6, 700, 8, 9}};

    array test = S*f;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],test[i])<<"i = "<<i ;
    }
}

TEST(Matrix,multiply_linewise){
    const Matrix S(RelaxationPar(1,10,100));
    const array f = {{1,2,3,4,5,6,7,8,9}};
    // const array vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};
    const array vergleich = {{ 1, 2, 30, 4, 500, 6, 700, 8, 9}};

    array test;

    for(int i = 0; i<9;i++)
    {
        test[i] = S.linewise(f,i);
        EXPECT_DOUBLE_EQ(vergleich[i],test[i])<<"i = "<<i ;
    }
}

TEST(Matrix,identity){
    const array f = {{1,2,3,4,5,6,7,8,9}};
    const array f0 = {{0,0,0,0,0,0,0,0,0}};

    const Matrix Identity = Matrix(true);
    const Matrix Zeros = Matrix();
    array test;

    EXPECT_EQ(f, Identity * f);
    EXPECT_EQ(f0, Zeros * f);
    for(int i = 0; i<9;i++)
    {
        test[i] = Identity.linewise(f,i);
    }
    EXPECT_EQ(f, test);
}

TEST(Matrix,plus_times){
    const Matrix Identity = Matrix(true);
    
    EXPECT_EQ(Identity+Identity, Identity*2);

    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX+TRAFO_MATRIX, TRAFO_MATRIX*3);
    EXPECT_EQ(TRAFO_MATRIX+TRAFO_MATRIX, (TRAFO_MATRIX*3)-TRAFO_MATRIX);
}

TEST(MRT,trafo){
    /// testet ob die Differenz im Geschw.-Raum gleich der Rücktransformierten Differenz im moment-Raum ist
    ParamSet param;
    DistributionSetType phi = param.getPhi();
    const array f = {{1,2,3,4,5,6,7,8,9}};
    Cell testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();

    VeloSet u = testCell.getU();
    DistributionSetType fEq  = eqDistro(rho,u,phi);
    DistributionSetType vergleich;
    vergleich[0] = array_diff(f, fEq[0]);
    vergleich[1] = array_diff(f, fEq[1]);
    array m = TRAFO_MATRIX * f;
    DistributionSetType mEq;
    mEq[0] = TRAFO_MATRIX * fEq[0];
    mEq[1] = TRAFO_MATRIX * fEq[1];
    DistributionSetType transformed;
    transformed[0] = INV_TRAFO_MATRIX * array_diff(m,mEq[0]);
    transformed[1] = INV_TRAFO_MATRIX * array_diff(m,mEq[1]);

    for(int i = 0; i<9;i++)
    {
        EXPECT_NEAR(vergleich[0][i],transformed[0][i],1e-10) ;
        EXPECT_NEAR(vergleich[1][i],transformed[1][i],1e-10) ;
    }
}

TEST(MRT,mass){
    ParamSet param;
    DistributionSetType phi = param.getPhi();
    const array f = {{1,2,3,4,5,6,7,8,9}};
    Cell testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();
    VeloSet u = testCell.getU();
    DistributionSetType fEq  = eqDistro(rho,u,phi);
    array m = TRAFO_MATRIX * f;
    array mEq = TRAFO_MATRIX * fEq[0];

    EXPECT_DOUBLE_EQ(m[0],mEq[0]);
    EXPECT_DOUBLE_EQ(m[3],mEq[3]);
    EXPECT_DOUBLE_EQ(m[5],mEq[5]);
}

TEST(ParamSet,Phi)
{
    ParamSet param(1,1,1,1000);
    DistributionSetType phi;
    phi = param.getPhi();

    array phiR = {{0.9992, 1.6e-4, 4e-5, 1.6e-4, 4e-5, 1.6e-4, 4e-5, 1.6e-4, 4e-5}};
    array phiB = {{0.2,0.16,0.04,0.16,0.04,0.16,0.04,0.16,0.04}};

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
    EXPECT_DOUBLE_EQ(50,newProcess.getReynoldsMax());
    EXPECT_DOUBLE_EQ(1e-3,newProcess.getMorton());
    EXPECT_DOUBLE_EQ(20,newProcess.getEotvos());
    EXPECT_DOUBLE_EQ(40,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1000,newProcess.getRhoL());
    EXPECT_DOUBLE_EQ(5,newProcess.getGamma());
    EXPECT_DOUBLE_EQ(0.1,newProcess.getDiameter());
    EXPECT_DOUBLE_EQ(10,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(1e-4,newProcess.getSigma());
    EXPECT_DOUBLE_EQ(10,newProcess.getGPhys());

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(0.6385640646055102,newProcess.getTau());
    EXPECT_DOUBLE_EQ(sqrt(3),newProcess.getSpeedlimit());
    EXPECT_DOUBLE_EQ(0.0025,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(0.00014433756729740645,newProcess.getTimestep());
    EXPECT_DOUBLE_EQ(0.002,newProcess.getNu());
    EXPECT_DOUBLE_EQ(800,newProcess.getDelRho());
}

TEST(Preprocess,FileInput){
    Preprocess newProcess = read_preprocess_file("preprocessFile");

    // test the given parameters (default values)
    EXPECT_DOUBLE_EQ(75,newProcess.getReynoldsMax());
    EXPECT_DOUBLE_EQ(0.0015,newProcess.getMorton());
    EXPECT_DOUBLE_EQ(25,newProcess.getEotvos());
    EXPECT_DOUBLE_EQ(45,newProcess.getResolution());
    EXPECT_DOUBLE_EQ(1200,newProcess.getRhoL());
    EXPECT_DOUBLE_EQ(2.75,newProcess.getGamma());
    EXPECT_DOUBLE_EQ(0.1,newProcess.getDiameter());
    EXPECT_DOUBLE_EQ(7.177827488211825,newProcess.getSoundspeed());
    EXPECT_DOUBLE_EQ(0.00291447961305505,newProcess.getSigma());
    EXPECT_DOUBLE_EQ(8.967172490810192,newProcess.getGPhys());   
 }

// TEST(Preprocess,refine){
//     Preprocess newProcess;
//     const double Reynolds = newProcess.getReynoldsMax();
//     const double c_s = newProcess.getSoundspeed();         
//     const double tau_r = newProcess.getTau()-0.5;
//     const double speedlimit = newProcess.getSpeedlimit();     
//     const double spacestep = newProcess.getSpacestep();  
//     const double timestep = newProcess.getTimestep();   
//     const double delRho = newProcess.getDelRho();     
//     const double nu = newProcess.getNu();

//     newProcess.refine();

//     EXPECT_DOUBLE_EQ(Reynolds*1.1, newProcess.getReynoldsMax());
//     EXPECT_DOUBLE_EQ(c_s*1.1, newProcess.getSoundspeed());
//     EXPECT_DOUBLE_EQ(tau_r/1.1 ,newProcess.getTau()-0.5);
//     EXPECT_DOUBLE_EQ(timestep / 1.1, newProcess.getTimestep());
//     EXPECT_DOUBLE_EQ(speedlimit*1.1, newProcess.getSpeedlimit());
//     EXPECT_DOUBLE_EQ(spacestep,newProcess.getSpacestep());
//     EXPECT_DOUBLE_EQ(delRho, newProcess.getDelRho());
//     EXPECT_DOUBLE_EQ(nu,newProcess.getNu());
// }

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
    EXPECT_DOUBLE_EQ(0.00084327404271156506,params.getSigma());
    EXPECT_DOUBLE_EQ(0.00014491665424627533,params.getG());
     }

TEST(timetrack,basic){
    Timetrack track(0.605,1.1);
    EXPECT_DOUBLE_EQ(0, track.getTime());
    track.timestep();
    track.timestep();
    track.timestep();
    EXPECT_DOUBLE_EQ(1.815, track.getTime());
    track.refine();
    EXPECT_DOUBLE_EQ(1.815, track.getTime());
    track.timestep();
    track.timestep();
    EXPECT_DOUBLE_EQ(2.915, track.getTime());
    track.refine();
    track.timestep();
    track.timestep();
    track.timestep();
    EXPECT_DOUBLE_EQ(4.415, track.getTime());
}

TEST(timetrack,FileInput){
    Preprocess newProcess = read_preprocess_file("preprocessFile");
    Timetrack newTimetrack = read_timetrack_file(newProcess, "preprocessFile");
     // test the given parameters 
    EXPECT_DOUBLE_EQ(newProcess.getTimestep(),newTimetrack.getDTini());
    EXPECT_DOUBLE_EQ(1.15,newTimetrack.getFactor());
    EXPECT_DOUBLE_EQ(4e5,newTimetrack.getMaxCount());

    EXPECT_DOUBLE_EQ(100,newTimetrack.getTechPlotInt());
    EXPECT_DOUBLE_EQ(1000,newTimetrack.getRestartInt());
}

TEST(Vector,scalar){
    Vector v0, v1(1,2), v2(3,4);
    EXPECT_DOUBLE_EQ(0, v0*v1);
    EXPECT_DOUBLE_EQ(0, v2*v0);
    EXPECT_DOUBLE_EQ(11, v1*v2);
    EXPECT_DOUBLE_EQ(11, v2*v1);
}

TEST(Vector,angle){
    Vector g(1,1);

    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION[0]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION[1]));
    EXPECT_DOUBLE_EQ(1, g.Angle(DIRECTION[2]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.Angle(DIRECTION[3]));
    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION[4]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION[5]));
    EXPECT_DOUBLE_EQ(-1, g.Angle(DIRECTION[6]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.Angle(DIRECTION[7]));
    EXPECT_DOUBLE_EQ(0, g.Angle(DIRECTION[8]));

    Vector g1(1e-10, -1e-10), g2(1e-6, -1e-6);
    EXPECT_DOUBLE_EQ(1,g1.Angle(g2));
    EXPECT_DOUBLE_EQ(0,g1.Angle(DIRECTION[0]));
}

#endif