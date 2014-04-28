#include"Lattice.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice::Lattice(int x_size, int y_size,double fzero_dense, double fzero_dilute):

xsize(x_size)
,ysize(y_size)
,data(new field(boost::extents[xsize][ysize]))
,param()
{
    for (int x = 0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y] = Cell(fzero_dense,fzero_dilute);
        }
    }
}

Lattice::Lattice(const Lattice& other):
xsize(other.getSize()[0])
,ysize(other.getSize()[1])
,data(new field(boost::extents[xsize][ysize]))
,param(other.getParams())
{
    (*data) = other.getData();
}


Lattice::~Lattice(){
    delete data;
    data = NULL;
}

//=========================== OPERATIONS ===========================

void Lattice::equilibriumIni()
{
    Cell tmp;
    DistributionSetType eqDis;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = (*data)[i][j];
            tmp.calcRho();
            ColSet rho = tmp.getRho();
            VeloSet u = tmp.getU();
            eqDis = eqDistro(rho,u,param.getPhi());
            tmp.setF(eqDis);
            (*data)[i][j] = tmp;
        }
    }
}

void Lattice::balance(double& mass, double& momentum)const
{
    VeloSet u;
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
            momentum += rho * sqrt(u[0]*u[0]);
        }
    }
}

void Lattice::mass_balance(double& liquid_mass, double& gas_mass)const
{
    ColSet rho;
    liquid_mass = 0;
    gas_mass = 0;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
            rho = (*data)[i][j].getRho();
            liquid_mass += rho[0];
            gas_mass += rho[1];
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
        tmp = x + DIRECTION[q].x;
        if (tmp<0) tmp += xsize;
        if (tmp>= xsize) tmp -= xsize;
        dir[q].x = tmp;

        tmp = y + DIRECTION[q].y;
        if (tmp<0) tmp += ysize;
        if (tmp>= ysize) tmp -= ysize;
        dir[q].y = tmp;
    }
    return dir;
}

const Vector Lattice::getGradient(int x, int y)const
{
    Vector grad(0,0);
    double tmpDelta;

    const direction dir = directions(x,y);
    for (int q=0;q<13;q++)
    {
        tmpDelta = GRAD_WEIGHTS[q] * (*data)[ dir[q].x ][ dir[q].y ].getDeltaRho();
        grad.x += DIRECTION[q].x * tmpDelta;
        grad.y += DIRECTION[q].y * tmpDelta;
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

bool Lattice::collideAll(int threads, bool gravity, bool isLimitActive)
{
    bool success(true);
    field *newData = new field(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);

    const double beta = param.getBeta();
    const DistributionSetType phi = param.getPhi();
    const int range = xsize * ysize;
    const double rhoRedFixed = param.getRhoR();
    const RelaxationPar relax = param.getRelaxation();
    const double dt = param.getDeltaT();
    const double speedlimit = param.getSpeedlimit();

    double g(0);
    if(gravity == true) g = param.getG();
    // if(gravity == true) g = 1e-1;


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
                DistributionSetType  fTmp;
                const DistributionSetType fCell = tmpCell.getF();

                const ColSet rho_k = tmpCell.getRho();
                const double rho = sum(rho_k);

                const Vector G(0 ,  g*(rho - rho_k[0]));
                // const Vector G(0 , - rho * g);

                VeloSet u = tmpCell.getU();

                if (isLimitActive == true)
                {
                    // check for exceptions
                    if ( u[0].Abs() > speedlimit || u[1].Abs() > speedlimit)
                    {
                        #pragma omp critical(Output)
                        cout<<"u[0].Abs() = "<< u[0].Abs() << "\tu[1].Abs() = "<< u[1].Abs() << "\tspeedlimit = " << speedlimit << endl;

                        #pragma omp critical(ExceptionLike)
                        success = false;
                    } 
                }   // end if(isLimitActive == true) 

                if(gravity == true){
                u[0] = u[0] + G *  (dt/(2* rho)) ;
                u[1] = u[1] + G *  (dt/(2* rho)) ;
                }

                const DistributionSetType fEq = eqDistro(rho_k, u, phi);
                const DistributionSetType diff = distro_diff(fCell, fEq);
            
                const double omega = param.getOmega(tmpCell.calcPsi());
                const Matrix relaxation_matrix(relax,omega);
                
                const Matrix forcing_factor = Matrix(true) - (relaxation_matrix*0.5);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
                const DistributionSetType first_forcing_term = forcing_factor * (TRAFO_MATRIX * calculate_forcing_term(G,u)); // F' = (I - 0.5 S) M F
                const DistributionSetType second_forcing_term = INV_TRAFO_MATRIX * first_forcing_term;    // M^{-1} F'

                const DistributionSetType single_phase_col = INV_TRAFO_MATRIX * (relaxation_matrix * (TRAFO_MATRIX * diff));
                          
                const ColSet A_k = param.getAk(omega);
                const Vector grad = getGradient(x,y);
                const double av = grad.Abs();

                double scal;
                double fges;
                double recolor;
                double final_forcing_term(0);                

                for (int q=0; q<9; q++)
                {
                    // gradient based two phase
                    scal = grad*DIRECTION[q];

                    double gradient_collision(0);
                    if (av > 0 ) gradient_collision = av/2 * (WEIGHTS[q] * ( scal*scal )/(av*av) - B[q]);

                    for (int color=0;color<=1; color++)
                    {
                        // if (gravity == true) final_forcing_term = second_forcing_term[color][q] ;
                        // fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q] + dt * final_forcing_term + A_k[color] * gradient_collision;
                        fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];// + dt * final_forcing_term;
                        // if (gravity == true) fTmp[color][q] +=  dt * second_forcing_term[color][q];
                    }
                    
                    for (int color=0;color<=1; color++)
                    {
                        fTmp[color][q] += A_k[color] * gradient_collision;
                    }

                    fges = fTmp[0][q]+fTmp[1][q];

                    // recoloring
                    if (rho > 0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) *  grad.Angle(DIRECTION[q])   * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q)); 
                    else recolor = 0;

                    if (rho > 0)
                    {
                        fTmp[0][q] = rho_k[0]/rho  * fges + recolor;
                        fTmp[1][q] = rho_k[1]/rho * fges - recolor;
                    }

                    if (gravity == true){
                        fTmp[0][q] +=  dt * second_forcing_term[0][q];
                        fTmp[1][q] +=  dt * second_forcing_term[1][q];
                    } 

                    // if (fTmp[0][q] < 0) fTmp[0][q] = 0;
                    // if (fTmp[1][q] < 0) fTmp[1][q] = 0;

                } // end for 

                tmpCell.setF(fTmp);
                tmpCell.calcRho();
            }
            #pragma omp critical(Zuweisung2)
            (*newData)[x][y] = tmpCell;
        }
    }
    if(success == true){
        delete data;
        data = newData;
    }
    else{
        delete newData;
    }     

    return success;
}

void Lattice::closedBox()
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

void Lattice::bottomWall()
{
    const Cell wall(0,0,true);

    for (int x=0; x<xsize; x++)
    {
        (*data)[x][0] = wall;
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

//=========================== ACCESSORS ===========================

const ColSet Lattice::getSize()const
{
    ColSet pony = {{xsize, ysize}};
    return pony;
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
    DistributionSetType f = (*data)[x][y].getF();
    f[color] = nf;
    (*data)[x][y].setF(f);
}

void Lattice::setF(int x, int y, int color, int pos, double value)
{
    DistributionSetType f = (*data)[x][y].getF();
    f[color][pos] = value;
    (*data)[x][y].setF(f);
}

//=========================== OPERATOR ===========================

Lattice& Lattice::operator=(const Lattice& other){
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1]);
    this->setParams(other.getParams());

    return *this;
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

///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

inline void Lattice::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}

void Lattice::streamAndBouncePull(Cell& tCell, const direction& dir)const
{
    const DistributionSetType f = tCell.getF();
    DistributionSetType ftmp;
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

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const DistributionSetType eqDistro(const ColSet& rho_k, const VeloSet& u, const DistributionSetType& phi)
{
    DistributionSetType feq;
    Vector velo = (u[0] * rho_k[0] + u[1] * rho_k[1]) * (1.0 / (rho_k[0]+rho_k[1])) ;
    double usqr = velo*velo;
   
    for (int i=0; i<9; i++)
    {
        double scal = velo*DIRECTION[i];
        for (int color = 0; color<=1; color++)
        {
            feq[color][i] = rho_k[color] * ( phi[color][i] + WEIGHTS[i] * ( 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
        }
    }
    return feq;
}

const DistributionSetType calculate_forcing_term(Vector G, VeloSet u)
{
    DistributionSetType forcing_term;
    for (int i=0; i<9; i++)
    {
        forcing_term[0][i] = 0;
        forcing_term[1][i] = WEIGHTS[i] * (G * ( DIRECTION[i] * (DIRECTION[i] * u[1]) + DIRECTION[i] - u[1])) ;
    }
    return forcing_term;
}