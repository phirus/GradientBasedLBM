#include"Lattice2D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice2D::Lattice2D(int x_size, int y_size,double fzero_dense, double fzero_dilute):
xsize(x_size)
,ysize(y_size)
,data(new field2D(boost::extents[xsize][ysize]))
,param()
{
    for (int x = 0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y] = Cell2D(fzero_dense,fzero_dilute);
        }
    }
}

Lattice2D::Lattice2D(const Lattice2D& other):
xsize(other.getSize()[0])
,ysize(other.getSize()[1])
,data(new field2D(boost::extents[xsize][ysize]))
,param(other.getParams())
{
    (*data) = other.getData();
}


Lattice2D::~Lattice2D(){
    delete data;
    data = NULL;
}

//=========================== OPERATIONS ===========================

void Lattice2D::equilibriumIni()
{
    Cell2D tmp;
    DistributionSetType2D eqDis;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = (*data)[i][j];
            tmp.calcRho();
            ColSet rho = tmp.getRho();
            VeloSet2D u = tmp.getU();
            eqDis = eqDistro(rho,u,param.getPhi2D());
            tmp.setF(eqDis);
            (*data)[i][j] = tmp;
        }
    }
}

void Lattice2D::balance(double& mass, double& momentum)const
{
    VeloSet2D u;
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

void Lattice2D::mass_balance(double& liquid_mass, double& gas_mass)const
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

void Lattice2D::overallRho()
{
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
        }
    }
}

direction2D Lattice2D::directions(int x, int y)const
{
    direction2D dir;
    int tmp;
    for (int q=0; q<13; q++)
    {
        tmp = x + DIRECTION_2D[q].x;
        if (tmp<0) tmp += xsize;
        if (tmp>= xsize) tmp -= xsize;
        dir[q].x = tmp;

        tmp = y + DIRECTION_2D[q].y;
        if (tmp<0) tmp += ysize;
        if (tmp>= ysize) tmp -= ysize;
        dir[q].y = tmp;
    }
    return dir;
}

const Vector2D Lattice2D::getGradient(int x, int y)const
{
    Vector2D grad(0,0);
    double tmpDelta;

    const direction2D dir = directions(x,y);
    for (int q=0;q<13;q++)
    {
        tmpDelta = GRAD_WEIGHTS_2D[q] * (*data)[ dir[q].x ][ dir[q].y ].getDeltaRho();
        grad.x += DIRECTION_2D[q].x * tmpDelta;
        grad.y += DIRECTION_2D[q].y * tmpDelta;
    }
    return grad;
}

void Lattice2D::streamAll(int threads)
{
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);
    const int range = xsize * ysize;

    #pragma omp parallel
    {
        #pragma omp for        
        for (int index = 0;  index < range; index++)
        {
            int x,y;
            linearIndex(index,x,y);
            const direction2D dir = directions(x,y);
            Cell2D tmpCell = (*data)[x][y];

            if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

            tmpCell.calcRho();
            #pragma omp critical(Zuweisung)
            (*newData)[x][y] = tmpCell;
        }
    }
    delete data;
    data = newData;
}

bool Lattice2D::collideAll(int threads, bool gravity, bool isLimitActive)
{
    bool success(true);
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);

    const double beta = param.getBeta();
    const DistributionSetType2D phi = param.getPhi2D();
    const int range = xsize * ysize;
    // const double rhoRedFixed = param.getRhoR();
    const RelaxationPar2D relax = param.getRelaxation2D();
    const double dt = param.getDeltaT();
    // const double speedlimit = param.getSpeedlimit();

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
            Cell2D tmpCell = (*data)[x][y];

            if (tmpCell.getIsSolid() == false)
            {
                DistributionSetType2D  fTmp;
                const DistributionSetType2D fCell = tmpCell.getF();

                const ColSet rho_k = tmpCell.getRho();
                const double rho = sum(rho_k);

                const Vector2D G(0 ,  g*(rho - rho_k[0]));
                // const Vector2D G(0 , - rho * g);

                VeloSet2D u = tmpCell.getU();

                // if (isLimitActive == true)
                // {
                //     // check for exceptions
                //     if ( u[0].Abs() > speedlimit || u[1].Abs() > speedlimit)
                //     {
                //         #pragma omp critical(Output)
                //         cout<<"u[0].Abs() = "<< u[0].Abs() << "\tu[1].Abs() = "<< u[1].Abs() << "\tspeedlimit = " << speedlimit << endl;

                //         #pragma omp critical(ExceptionLike)
                //         success = false;
                //     } 
                // }   // end if(isLimitActive == true) 

                if(gravity == true){
                u[0] = u[0] + G *  (dt/(2* rho)) ;
                u[1] = u[1] + G *  (dt/(2* rho)) ;
                }

                const DistributionSetType2D fEq = eqDistro(rho_k, u, phi);
                const DistributionSetType2D diff = distro_diff_2D(fCell, fEq);
            
                const double omega = param.getOmega(tmpCell.calcPsi());
                const Matrix2D relaxation_matrix(relax,omega);
                
                const Matrix2D forcing_factor = Matrix2D(true) - (relaxation_matrix*0.5);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
                const DistributionSetType2D first_forcing_term = forcing_factor * (TRAFO_MATRIX2D * calculate_forcing_term(G,u)); // F' = (I - 0.5 S) M F
                const DistributionSetType2D second_forcing_term = INV_TRAFO_MATRIX2D * first_forcing_term;    // M^{-1} F'

                const DistributionSetType2D single_phase_col = INV_TRAFO_MATRIX2D * (relaxation_matrix * (TRAFO_MATRIX2D * diff));
                          
                const ColSet A_k = param.getAk(omega);
                const Vector2D grad = getGradient(x,y);
                const double av = grad.Abs();

                double scal;
                double fges;
                double recolor;
                // double final_forcing_term(0);                

                for (int q=0; q<9; q++)
                {
                    // gradient based two phase
                    scal = grad*DIRECTION_2D[q];

                    double gradient_collision(0);
                    if (av > 0 ) gradient_collision = av/2 * (WEIGHTS_2D[q] * ( scal*scal )/(av*av) - B_2D[q]);

                    for (int color=0;color<=1; color++)
                    {
                        // if (gravity == true) final_forcing_term = second_forcing_term[color][q] ;
                        // fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q] + dt * final_forcing_term + A_k[color] * gradient_collision;
                        fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];// + dt * final_forcing_term;
                        if (gravity == true) fTmp[color][q] +=  dt * second_forcing_term[color][q];
                    }
                    
                    for (int color=0;color<=1; color++)
                    {
                        fTmp[color][q] += A_k[color] * gradient_collision;
                    }

                    fges = fTmp[0][q]+fTmp[1][q];

                    // recoloring
                    if (rho > 0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) *  grad.Angle(DIRECTION_2D[q])   * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q)); 
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

void Lattice2D::closedBox()
{
    const Cell2D wall(0,0,true);

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

void Lattice2D::bottomWall()
{
    const Cell2D wall(0,0,true);

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

void Lattice2D::genericWall(std::vector<double> x, std::vector<double> y, const Vector2D& u_w)
{
    if(x.size() == y.size())
    {
        Cell2D wall(0,0,true);
        wall.setSolidVelocity(u_w);
        
        for(unsigned int i = 0; i< x.size(); i++)
        {
            (*data)[x[i]][y[i]] = wall;
        }
    }
    else throw "vector size mismatch";    
}

void Lattice2D::lidDrivenCavity(const Vector2D& u_w)
{
    Cell2D wall(0,0,true);

    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;   // left wall
        (*data)[xsize-1][y] = wall; // right wall
    }

    for (int x=0; x<xsize; x++) 
    {
        (*data)[x][0] = wall; // bottom wall
    }

    wall.setSolidVelocity(u_w);
    for (int x=0; x<xsize; x++) 
    {
         (*data)[x][ysize-1] = wall;    // top wall
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][y].calcRho();
        }
    }
}

void Lattice2D::shearWall(const Vector2D& u_w)
{
    Cell2D wall(0,0,true);

    for (int y=0; y<ysize; y++)
    {
        (*data)[xsize-1][y] = wall; // right wall
    }

    wall.setSolidVelocity(u_w);
    for (int y=0; y<ysize; y++)
    {
        (*data)[0][y] = wall;   // left wall
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

const ColSet Lattice2D::getSize()const
{
    ColSet pony = {{double(xsize), double(ysize)}}; // need to cast int -> double
    return pony;
}

void Lattice2D::setData(const field2D& ndata, int x, int y){
    data->resize(boost::extents[x][y]);
    *data = ndata;
    xsize = x;
    ysize = y;
}

void Lattice2D::setCell(int x, int y, const Cell2D& ncell)
{
    if (y >= 0 && y < ysize && x >= 0 && x < xsize) (*data)[x][y] = ncell;
}

void Lattice2D::setF(int x, int y, int color, const array2D& nf)
{
    DistributionSetType2D f = (*data)[x][y].getF();
    f[color] = nf;
    (*data)[x][y].setF(f);
}

void Lattice2D::setF(int x, int y, int color, int pos, double value)
{
    DistributionSetType2D f = (*data)[x][y].getF();
    f[color][pos] = value;
    (*data)[x][y].setF(f);
}

//=========================== OPERATOR ===========================

Lattice2D& Lattice2D::operator=(const Lattice2D& other){
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1]);
    this->setParams(other.getParams());

    return *this;
}

const bool Lattice2D::operator==(const Lattice2D& other)const
{
    bool exit = true;
    ColSet extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        field2D otherData = other.getData();
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

inline void Lattice2D::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}

void Lattice2D::streamAndBouncePull(Cell2D& tCell, const direction2D& dir)const
{
    const DistributionSetType2D f = tCell.getF();
    DistributionSetType2D ftmp;
    for (int color = 0; color<=1;color++)
    {
        ftmp[color][0] = (*data)[ dir[0].x ][ dir[0].y ].getF()[color][0];

        for (int i=1;i<9;i++)
        {
            const Cell2D neighbor = (*data)[ dir[PULL_INDEX_2D[i]].x ][ dir[PULL_INDEX_2D[i]].y ];
            if (neighbor.getIsSolid() == false) // if(neighbor not solid?) -> stream
            {
                ftmp[color][i] = neighbor.getF()[color][i];    
            }
            else  // else -> bounce back
            {
                const ColSet rho = tCell.getRho();
                ftmp[color][i] = f[color][PULL_INDEX_2D[i]] - (2.0 * WEIGHTS_2D[i] * rho[color] * (DIRECTION_2D[i] * neighbor.getU()[0]) ) ;
            } 
        } // end for i
    } // end for color
    tCell.setF(ftmp);
}

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const DistributionSetType2D eqDistro(const ColSet& rho_k, const VeloSet2D& u, const DistributionSetType2D& phi)
{
    DistributionSetType2D feq;
    Vector2D velo = (u[0] * rho_k[0] + u[1] * rho_k[1]) * (1.0 / (rho_k[0]+rho_k[1])) ;
    double usqr = velo*velo;
   
    for (int i=0; i<9; i++)
    {
        double scal = velo*DIRECTION_2D[i];
        for (int color = 0; color<=1; color++)
        {
            feq[color][i] = rho_k[color] * ( phi[color][i] + WEIGHTS_2D[i] * ( 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
        }
    }
    return feq;
}

const DistributionSetType2D calculate_forcing_term(Vector2D G, VeloSet2D u)
{
    DistributionSetType2D forcing_term;
    for (int i=0; i<9; i++)
    {
        forcing_term[0][i] = 0;
        forcing_term[1][i] = WEIGHTS_2D[i] * (G * ( DIRECTION_2D[i] * (DIRECTION_2D[i] * u[1]) + DIRECTION_2D[i] - u[1])) ;
    }
    return forcing_term;
}