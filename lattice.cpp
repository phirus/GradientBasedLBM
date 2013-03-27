#include "lattice.h"

Lattice::Lattice(int x_size, int y_size,double fzero_red, double fzero_blue):data(new field(boost::extents[x_size][y_size])),param()
{
    xsize = x_size;
    ysize = y_size;
    for (int x = 0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y] = Cell(fzero_red,fzero_blue);
        }
    }
}
Lattice::~Lattice(){
    delete data;
    data = NULL;
}

void Lattice::setData(const field& ndata, int x, int y){
    data->resize(boost::extents[x][y]);
    *data = ndata;
    xsize = x;
    ysize = y;
}

void Lattice::setCell(int x, int y, const Cell& ncell)
{
    if (y >= 0 && y < ysize && x >= 0 && x < xsize) (*data)[x][y] = ncell;
}

void Lattice::setF(int x, int y, int color, const array& nf)
{
    FSet f = (*data)[x][y].getF();
    f[color] = nf;
    (*data)[x][y].setF(f);
}

void Lattice::setF(int x, int y, int color, int pos, double value)
{
    FSet f = (*data)[x][y].getF();
    f[color][pos] = value;
    (*data)[x][y].setF(f);
}

const ColSet Lattice::getSize()const
{
    ColSet pony = {{xsize, ysize}};
    return pony;
}

void Lattice::initialize()
{
    const Cell wall(0,0,true);

    for (int x=0; x<xsize; x++)
    {
        (*data)[x][0] = wall;
        (*data)[x][ysize-1] = wall;
    }
    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;
        (*data)[xsize-1][y] = wall;
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice::equilibriumIni()
{
    Cell tmp;
    FSet eqDis;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = (*data)[i][j];
            tmp.calcRho();
            ColSet rho = tmp.getRho();
            Vector u = tmp.getU();
            eqDis = eqDistro(rho,u,param.getPhi());
            tmp.setF(eqDis);
            (*data)[i][j] = tmp;
        }
    }
}

void Lattice::balance(double& mass, double& momentum)const
{
    Vector u;
    double rho;

    mass = 0;
    momentum = 0;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
            rho = sum((*data)[i][j].getRho());
            u = (*data)[i][j].getU();

            mass += rho;
            momentum += rho * sqrt(u*u);
        }
    }
}

void Lattice::overallRho()
{
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
        }
    }
}

direction Lattice::directions(int x, int y)const
{
    direction dir;
    int tmp;
    for (int q=0; q<13; q++)
    {
        tmp = x + e[q].x;
        if (tmp<0) tmp += xsize;
        if (tmp>= xsize) tmp -= xsize;
        dir[q].x = tmp;

        tmp = y + e[q].y;
        if (tmp<0) tmp += ysize;
        if (tmp>= ysize) tmp -= ysize;
        dir[q].y = tmp;
    }
    return dir;
}

const Vector Lattice::getGradient(int x, int y)const
{
    Vector grad(0,0);;
    double tmpDelta;

    const direction dir = directions(x,y);
    for (int q=0;q<13;q++)
    {
        tmpDelta = xi[q] * (*data)[ dir[q].x ][ dir[q].y ].getDeltaRho();
        grad.x += e[q].x * tmpDelta;
        grad.y += e[q].y * tmpDelta;
    }
    return grad;
}

void Lattice::streamAll(int threads)
{
    field *newData = new field(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);
    const int range = xsize * ysize;
#pragma omp parallel
    {
#pragma omp for
        for (int index = 0;  index < range; index++)
        {
            int x,y;
            linearIndex(index,x,y);
            const direction dir = directions(x,y);
            Cell tmpCell = (*data)[x][y];

            if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

            tmpCell.calcRho();
#pragma omp critical(Zuweisung)
            (*newData)[x][y] = tmpCell;
        }
    }
    delete data;
    data = newData;
}

void Lattice::collideAll(int threads, bool gravity)
{
    omp_set_num_threads (threads);

    field *newData = new field(boost::extents[xsize][ysize]);

    const double beta = param.getBeta();
    const FSet phi = param.getPhi();
    const double g = param.getG();
    const int range = xsize * ysize;
    const double rhoRedFixed = param.getRhoR();
    const RelaxationPar relax = param.getRelaxation();
    const double dt = param.getDeltaT();

#pragma omp parallel
    {
#pragma omp for
        for (int index = 0;  index < range; index++)
        {
            int x,y;
            linearIndex(index,x,y);
            Cell tmpCell = (*data)[x][y];

            if (tmpCell.getIsSolid() == false)
            {
                FSet  fTmp;
                FSet fCell = tmpCell.getF();
                FSet Diff;

                const ColSet rho_k = tmpCell.getRho();
                const double rho = sum(rho_k);

//                const Vector G(0 ,  g*(rhoRedFixed - rho));
                const Vector G(0 , - rho * g);

//                const Vector u = tmpCell.calcU();
                const Vector u = tmpCell.getU() + G *  (dt/(2* rho)) ;

                const FSet fEq = eqDistro(rho_k, u, phi);

                Diff[0] = arrayDiff(fCell[0],fEq[0]);
                Diff[1] = arrayDiff(fCell[1],fEq[1]);

                const double omega = param.getOmega(tmpCell.calcPsi());

                const ColSet A_k = param.getAk(omega);

                const Vector grad = getGradient(x,y);
                const double av = grad.abs();

                double two_phase;
                double scal;
                double fges;
                double recolor;
                double forcingTerm = 0;

                for (int q=0; q<9; q++)
                {
                    scal = grad*e[q];
                    if (av > 0) two_phase = av/2 * (w[q] * ( scal*scal )/(av*av) - B[q]);
                    else two_phase = 0;

                    if (gravity == true) forcingTerm = (1- 0.5*omega) * w[q] * (G * ( e[q] * (e[q] * u) + e[q] - u )) ;

                    for (int color=0;color<=1; color++)
                    {
                        fTmp[color][q] =  fCell[color][q] - omega * Diff[color][q] + A_k[color] * two_phase + forcingTerm * dt;
                        if (fTmp[color][q] < 0) fTmp[color][q] = 0;
                    }
                    if (rho_k[0] > 0 && rho_k[1] > 0 && rho > 0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) *  grad.angle(e[q])   * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q));
                    else recolor = 0;

                    fges = fTmp[0][q]+fTmp[1][q];

                    if (rho > 0)
                    {
                        fTmp[0][q] = rho_k[0]/rho  * fges + recolor;
                        fTmp[1][q] = rho_k[1]/rho * fges - recolor;
                    }

                    if (fTmp[0][q] < 0) fTmp[0][q] = 0;
                    if (fTmp[1][q] < 0) fTmp[1][q] = 0;
                }
                tmpCell.setF(fTmp);
            }
#pragma omp critical(Zuweisung2)
            (*newData)[x][y] = tmpCell;
        }
    }
    delete data;
    data = newData;
}

const bool Lattice::operator==(const Lattice& other)const
{
    bool exit = true;
    ColSet extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        field otherData = other.getData();
        for (int x = 0; x< xsize;x++)
        {
            for (int y = 0; y<ysize;y++)
            {
                if (!((*data)[x][y]==otherData[x][y]))
                {
                    exit = false;
                    break;
                }
            }
        }
    }
    else exit = false;

    return exit;
}

inline void Lattice::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}

void Lattice::streamAndBouncePull(Cell& tCell, const direction& dir)const
{
    const FSet f = tCell.getF();
    FSet ftmp;
    for (int color = 0; color<=1;color++)
    {
        ftmp[color][0] = (*data)[ dir[0].x ][ dir[0].y ].getF()[color][0];

        if ((*data)[ dir[5].x ][ dir[5].y ].getIsSolid() == false)
        {
            ftmp[color][1] = (*data)[ dir[5].x ][ dir[5].y ].getF()[color][1];    // if(neighbor not solid?) -> stream
        }
        else ftmp[color][1] = f[color][5];                                                                                            // else -> bounce back

        if ((*data)[ dir[6].x ][ dir[6].y ].getIsSolid() == false)
        {
            ftmp[color][2] = (*data)[ dir[6].x ][ dir[6].y ].getF()[color][2];
        }
        else ftmp[color][2] = f[color][6];

        if ((*data)[ dir[7].x ][ dir[7].y ].getIsSolid() == false)
        {
            ftmp[color][3] = (*data)[ dir[7].x ][ dir[7].y ].getF()[color][3];
        }
        else ftmp[color][3] = f[color][7];

        if ((*data)[ dir[8].x ][ dir[8].y ].getIsSolid() == false)
        {
            ftmp[color][4] = (*data)[ dir[8].x ][ dir[8].y ].getF()[color][4];
        }
        else ftmp[color][4] = f[color][8];

        if ((*data)[ dir[1].x ][ dir[1].y ].getIsSolid() == false)
        {
            ftmp[color][5] = (*data)[ dir[1].x ][ dir[1].y ].getF()[color][5];
        }
        else ftmp[color][5] = f[color][1];

        if ((*data)[ dir[2].x ][ dir[2].y ].getIsSolid() == false)
        {
            ftmp[color][6] = (*data)[ dir[2].x ][ dir[2].y ].getF()[color][6];
        }
        else ftmp[color][6] = f[color][2];

        if ((*data)[ dir[3].x ][ dir[3].y ].getIsSolid() == false)
        {
            ftmp[color][7] = (*data)[ dir[3].x ][ dir[3].y ].getF()[color][7];
        }
        else ftmp[color][7] = f[color][3];

        if ((*data)[ dir[4].x ][ dir[4].y ].getIsSolid() == false)
        {
            ftmp[color][8] = (*data)[ dir[4].x ][ dir[4].y ].getF()[color][8];
        }
        else ftmp[color][8] = f[color][4];
    }
    tCell.setF(ftmp);
}

const FSet eqDistro(const ColSet& rho_k, const Vector& u, const FSet& phi)
{
    FSet feq;
    const double usqr = u*u;
    double scal;
    for (int i=0; i<9; i++)
    {
        scal = u*e[i];
        for (int color = 0; color<=1; color++)
        {
            feq[color][i] = rho_k[color] * ( phi[color][i] + w[i] * ( 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
        }
    }
    return feq;
}

const array arrayDiff(const array &one, const array &two)
{
    array a;
    for (int i=0; i<9; i++)
    {
        a[i] = one[i]-two[i];
    }
    return a;
}
