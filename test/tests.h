#ifndef TESTS_H
#define TESTS_H

#include"gtest/gtest.h"
#include"../src/lattice.h"
#include"../src/vector.h"
#include"../src/binaryIO.h"

using namespace std;

TEST(Cell,constructor0)
{
    array f = {{1,0,0,0,0,0,0,0,0}};
    Cell cell;
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_EQ(false,cell.getIsSolid());
    EXPECT_EQ(f,cell.getF()[1]);
}

TEST(Cell,constructorRed1)
{
    array f = {{1,0,0,0,0,0,0,0,0}};
    Cell cell(1,0);
    EXPECT_EQ(f,cell.getF()[0]);
    EXPECT_EQ(false,cell.getIsSolid());
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
    EXPECT_EQ(false,cell.getIsSolid());
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
    EXPECT_EQ(false,cell.getIsSolid());
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
    EXPECT_EQ(false,cell.getIsSolid());
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
    EXPECT_EQ(false,cell.getIsSolid());
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
    Vector u = cell.getU();
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
    FSet f = four.getF();
    f[1] = fr;
    four.setF(f);
    EXPECT_FALSE(one == four);
}

TEST(Matrix,trafo){
    const array verteilung = {{1,2,3,4,5,6,7,8,9}};
    const array vergleich = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array trafo = TrafoMatrix * verteilung;

    EXPECT_EQ(vergleich, trafo);
}

TEST(Matrix,backtrafo){
    const array vergleich = {{1,2,3,4,5,6,7,8,9}};
    const array verteilung = {{45,24,-12,-4,8,-12,0,-4,-4}};

    array backtrafo = InvTtrafoMatrix * verteilung;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i],backtrafo[i]);
    }
}

TEST(Matrix,relaxation_without_omega){
    const Matrix S(1,10,100);
    const array f = {{1,2,3,4,5,6,7,8,9}};
    const array vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};

    array test = S*f;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i]/3,test[i])<<"i = "<<i ;
    }

}

TEST(Matrix,relaxation_linewise){
    const Matrix S(1,10,100);
    const array f = {{1,2,3,4,5,6,7,8,9}};
    const array vergleich = {{ -48, -382, 194, 18, -206, 418, -206, 18, 194}};
    array test;

    for(int i = 0; i<9;i++)
    {
        test[i] = S.linewise(f,i);
        EXPECT_DOUBLE_EQ(vergleich[i]/3,test[i])<<"i = "<<i ;
    }
}

TEST(Matrix,relaxation_omega){
    Matrix S(1,10,20);
    S.addOmega(5);
    const array f = {{1,2,3,4,5,6,7,8,9}};
    const array vergleich = {{ -48, -77, 19, 33, -31, 83, -61, 33, 49}};

    array test = S*f;

    for(int i = 0; i<9;i++)
    {
        EXPECT_DOUBLE_EQ(vergleich[i]/3,test[i])<<"i = "<<i ;
    }

}

TEST(ParamSet,Phi)
{
    ParamSet param;
    FSet phi;
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
    EXPECT_FALSE(one == four);
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

    EXPECT_DOUBLE_EQ(0, g.angle(e[0]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.angle(e[1]));
    EXPECT_DOUBLE_EQ(1, g.angle(e[2]));
    EXPECT_DOUBLE_EQ(cos(PI/4), g.angle(e[3]));
    EXPECT_DOUBLE_EQ(0, g.angle(e[4]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.angle(e[5]));
    EXPECT_DOUBLE_EQ(-1, g.angle(e[6]));
    EXPECT_DOUBLE_EQ(-cos(PI/4), g.angle(e[7]));
    EXPECT_DOUBLE_EQ(0, g.angle(e[8]));

    Vector g1(1e-10, -1e-10), g2(1e-6, -1e-6);
    EXPECT_DOUBLE_EQ(1,g1.angle(g2));
    EXPECT_DOUBLE_EQ(0,g1.angle(e[0]));
}


TEST(Constants,BReis)
{
    Lattice lattice;

    EXPECT_DOUBLE_EQ(-4,B[0]*27);
    EXPECT_DOUBLE_EQ(2,B[1]*27);
    EXPECT_DOUBLE_EQ(5,B[2]*108);
}
TEST(Constants,W)
{
    Lattice lattice;

    EXPECT_DOUBLE_EQ(4,w.at(0)*9);
    EXPECT_DOUBLE_EQ(1,w.at(1)*9);
    EXPECT_DOUBLE_EQ(1,w.at(2)*36);
}

TEST(Constants,Xi)
{
    Lattice lattice;
    EXPECT_DOUBLE_EQ(0,xi.at(0));
    EXPECT_DOUBLE_EQ(32,xi.at(1)*120);
    EXPECT_DOUBLE_EQ(12,xi.at(2)*120);
    EXPECT_DOUBLE_EQ(1,xi.at(9)*120);
}

TEST(Lattice,constructor)
{
    Lattice lattice;
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
    Vector u = cell.getU();
    usqr = u*u;

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
    double abs = grad.abs();
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

TEST(MRT,trafo){
    /// testet ob die Differenz im Geschw.-Raum gleich der Rücktransformierten Differenz im moment-Raum ist
    ParamSet param;
    FSet phi = param.getPhi();
    const array f = {{1,2,3,4,5,6,7,8,9}};
    Cell testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();
    Vector u = testCell.getU();
    FSet fEq  = eqDistro(rho,u,phi);
    FSet vergleich;
    vergleich[0] = arrayDiff(f, fEq[0]);
    vergleich[1] = arrayDiff(f, fEq[1]);
    array m = TrafoMatrix * f;
    FSet mEq;
    mEq[0] = TrafoMatrix * fEq[0];
    mEq[1] = TrafoMatrix * fEq[1];
    FSet transformed;
    transformed[0] = InvTtrafoMatrix * arrayDiff(m,mEq[0]);
    transformed[1] = InvTtrafoMatrix * arrayDiff(m,mEq[1]);
    for(int i = 0; i<9;i++)
    {
        EXPECT_NEAR(vergleich[0][i],transformed[0][i],1e-10) ;
        EXPECT_NEAR(vergleich[1][i],transformed[1][i],1e-10) ;
    }
}
TEST(MRT,mass){
    ParamSet param;
    FSet phi = param.getPhi();
    const array f = {{1,2,3,4,5,6,7,8,9}};
    Cell testCell(f,f);
    testCell.calcRho();
    ColSet rho = testCell.getRho();
    Vector u = testCell.getU();
    FSet fEq  = eqDistro(rho,u,phi);
    array m = TrafoMatrix * f;
    array mEq = TrafoMatrix * fEq[0];

    EXPECT_DOUBLE_EQ(m[0],mEq[0]);
    EXPECT_DOUBLE_EQ(m[3],mEq[3]);
    EXPECT_DOUBLE_EQ(m[5],mEq[5]);
}

TEST(BinaryIO,output){
    Lattice lattice(150,150);
    binary_output(lattice);
    Lattice vergleich;
    EXPECT_FALSE(binary_input(vergleich,"existiertnicht.txt"));
    EXPECT_TRUE(binary_input(vergleich));
    EXPECT_EQ(lattice,vergleich);
}
 TEST(BinaryIO,paramLog){
     Lattice lattice(100,100);    
     EXPECT_NO_THROW(paramLogOut(lattice));
}
TEST(BinaryIO,queryTest){
    double value;
    EXPECT_FALSE( inputQuery("existiertnicht","test",value) );
    EXPECT_FALSE(inputQuery("queryTest","noflag",value));
    EXPECT_TRUE(inputQuery("queryTest","test",value));
    EXPECT_DOUBLE_EQ(13.4, value); 

}

TEST(BinaryIO,paramConfIn){
    ParamSet param(0.8, 1.4, 1.1, 1100);
    ParamSet input = getFileParams("paramInputTest");
    EXPECT_EQ(param, input);    

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

    // test the deduced parameters
    EXPECT_DOUBLE_EQ(0.6385640646055102,newProcess.getTau());
    EXPECT_DOUBLE_EQ(sqrt(3),newProcess.getSpeedlimit());
    EXPECT_DOUBLE_EQ(0.0025,newProcess.getSpacestep());
    EXPECT_DOUBLE_EQ(0.00014433756729740645,newProcess.getTimestep());
    EXPECT_DOUBLE_EQ(0.002,newProcess.getNu());
    EXPECT_DOUBLE_EQ(800,newProcess.getDelRho());
}



#endif
