#include"Lattice3D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice3D::Lattice3D(int x_size, int y_size, int z_size, double fzero_dense, double fzero_dilute):
xsize(x_size)
,ysize(y_size)
,zsize(z_size)
,data(new field3D(boost::extents[xsize][ysize][zsize]))
,param()
{
    for (int x = 0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            for (int z=0; z<zsize; z++)
            {
                (*data)[x][y][z] = Cell3D(fzero_dense,fzero_dilute);
            }
        }
    }
}

Lattice3D::Lattice3D(const Lattice3D& other):
xsize(other.getSize()[0])
,ysize(other.getSize()[1])
,zsize(other.getSize()[2])
,data(new field3D(boost::extents[xsize][ysize][zsize]))
,param(other.getParams())
{
    (*data) = other.getData();
}


Lattice3D::~Lattice3D(){
    delete data;
    data = NULL;
}

//=========================== OPERATIONS ===========================

void Lattice3D::equilibriumIni()
{
    Cell3D tmp;
    DistributionSetType3D eqDis;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            for (int k=0; k<zsize; k++)
            {
                tmp = (*data)[i][j][k];
                tmp.calcRho();
                ColSet rho = tmp.getRho();
                VeloSet3D u = tmp.getU();
                eqDis = eqDistro(rho,u,param.getPhi3D());
                tmp.setF(eqDis);
                (*data)[i][j][k] = tmp;
            }
        }
    }
}

void Lattice3D::balance(double& mass, double& momentum)const
{
    VeloSet3D u;
    double rho;

    mass = 0;
    momentum = 0;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            for (int k=0;k<zsize; k++)
            {
                (*data)[i][j][k].calcRho();
                rho = sum((*data)[i][j][k].getRho());
                u = (*data)[i][j][k].getU();

                mass += rho;
                momentum += rho * sqrt(u[0]*u[0]);
            }
        }
    }
}

void Lattice3D::mass_balance(double& liquid_mass, double& gas_mass)const
{
    ColSet rho;
    liquid_mass = 0;
    gas_mass = 0;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            for (int k=0; k<zsize; k++)
            {
                (*data)[i][j][k].calcRho();
                rho = (*data)[i][j][k].getRho();
                liquid_mass += rho[0];
                gas_mass += rho[1];
            }
        }
    }
}

void Lattice3D::overallRho()
{
    for (int i=0; i<xsize; i++)
    {
        for (int j=0; j<ysize; j++)
        {
            for (int k = 0; k<zsize; k++)
            {
                (*data)[i][j][k].calcRho();
            } // end for z
        } // end for y
    } // end for x
}

direction3D Lattice3D::directions(int x, int y, int z)const
{
    direction3D dir;
    int tmp;
    for (int q=0; q<33; q++)
    {
        tmp = x + DIRECTION_3D[q].x;
        if (tmp<0) tmp += xsize;
        if (tmp>= xsize) tmp -= xsize;
        dir[q].x = tmp;

        tmp = y + DIRECTION_3D[q].y;
        if (tmp<0) tmp += ysize;
        if (tmp>= ysize) tmp -= ysize;
        dir[q].y = tmp;

        tmp = z + DIRECTION_3D[q].z;
        if (tmp<0) tmp += zsize;
        if (tmp>= zsize) tmp -= zsize;
        dir[q].z = tmp;
    }
    return dir;
}

const Vector3D Lattice3D::getGradient(int x, int y, int z)const
{
    Vector3D grad(0,0,0);
    double tmpDelta;

    const direction3D dir = directions(x,y,z);
    for (int q=0;q<33;q++)
    {
        tmpDelta = GRAD_WEIGHTS_3D[q] * (*data)[ dir[q].x ][ dir[q].y ][ dir[q].z ].getDeltaRho();
        grad.x += DIRECTION_3D[q].x * tmpDelta;
        grad.y += DIRECTION_3D[q].y * tmpDelta;
        grad.z += DIRECTION_3D[q].z * tmpDelta;
    }
    return grad;
}

void Lattice3D::streamAll(int threads)
{
    field3D *newData = new field3D(boost::extents[xsize][ysize][zsize]);

    omp_set_num_threads (threads);
    const int range = xsize * ysize * zsize;

    #pragma omp parallel
    {
        #pragma omp for        
        for (int index = 0;  index < range; index++)
        {
            int x,y,z;
            linearIndex(index,x,y,z);
            const direction3D dir = directions(x,y,z);
            Cell3D tmpCell = (*data)[x][y][z];

            if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

            tmpCell.calcRho();
            #pragma omp critical(Zuweisung)
            (*newData)[x][y][z] = tmpCell;
        }
    }
    delete data;
    data = newData;
}

bool Lattice3D::collideAll(int threads, bool gravity, bool isLimitActive)
{
    bool success(true);
    field3D *newData = new field3D(boost::extents[xsize][ysize][zsize]);

     omp_set_num_threads (threads);

    const double beta = param.getBeta();
    const DistributionSetType3D phi = param.getPhi3D();
    const int range = xsize * ysize * zsize;
    // const double rhoRedFixed = param.getRhoR();
    const RelaxationPar3D relax = param.getRelaxation3D();
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
            int x,y,z;
            linearIndex(index,x,y,z);
            Cell3D tmpCell = (*data)[x][y][z];

            if (tmpCell.getIsSolid() == false)
            {
                DistributionSetType3D  fTmp;
                const DistributionSetType3D fCell = tmpCell.getF();

                const ColSet rho_k = tmpCell.getRho();
                const double rho = sum(rho_k);

                const Vector3D G(0 , 0, g*(rho - rho_k[0]));
                // const Vector2D G(0 , 0, - rho * g);

                VeloSet3D u = tmpCell.getU();

                if(gravity == true)
                {
                    u[0] = u[0] + G *  (dt/(2* rho)) ;
                    u[1] = u[1] + G *  (dt/(2* rho)) ;
                }

                const DistributionSetType3D fEq = eqDistro(rho_k, u, phi);
                const DistributionSetType3D diff = distro_diff_3D(fCell, fEq);
            
                const double omega = param.getOmega(tmpCell.calcPsi());
                const Matrix3D relaxation_matrix(relax,omega);
                
                const Matrix3D forcing_factor = Matrix3D(true) - (relaxation_matrix*0.5);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
                const DistributionSetType3D first_forcing_term = forcing_factor * (TRAFO_MATRIX3D * calculate_forcing_term(G,u)); // F' = (I - 0.5 S) M F
                const DistributionSetType3D second_forcing_term = INV_TRAFO_MATRIX3D * first_forcing_term;    // M^{-1} F'

                const DistributionSetType3D single_phase_col = INV_TRAFO_MATRIX3D * (relaxation_matrix * (TRAFO_MATRIX3D * diff));
                          
                const ColSet A_k = param.getAk(omega);
                const Vector3D grad = getGradient(x,y,z);
                const double av = grad.Abs();

                double scal;
                double fges;
                double recolor;
                // double final_forcing_term(0);                

                for (int q=0; q<19; q++)
                {
                    // gradient based two phase
                    scal = grad*DIRECTION_3D[q];

                    double gradient_collision(0);
                    if (av > 0 ) gradient_collision = av/2 * (WEIGHTS_3D[q] * ( scal*scal )/(av*av) - B_3D[q]);

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
                    if (rho > 0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) *  grad.Angle(DIRECTION_3D[q])   * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q)); 
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
            } // end if(solid)
            #pragma omp critical(Zuweisung2)
            (*newData)[x][y][z] = tmpCell;
        } // end for(index) 
    } // end #pragma parallel 
    if(success == true){
        delete data;
        data = newData;
    }
    else
    {
        delete newData;
    }  

    return success;
}

void Lattice3D::closedBox()
{
    const Cell3D wall(0,0,true);

    // top and bootom wall
    for (int x=0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y][0] = wall;
            (*data)[x][y][zsize-1] = wall;
        }
    }

    // left and right wall
    for (int z=0; z<zsize; z++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[0][y][z] = wall;
            (*data)[xsize-1][y][z] = wall;
        }
    }

    // back and front wall
    for (int z=0; z<zsize; z++)
    {
        for (int x=0; x<xsize; x++)
        {
            (*data)[x][0][z] = wall;
            (*data)[x][ysize-1][z] = wall;
        }
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            for (int z=0; z<zsize; z++)
            {
                (*data)[x][y][z].calcRho();
            }
        }
    }
}

void Lattice3D::bottomWall()
{
    const Cell3D wall(0,0,true);

    for (int x=0; x<xsize; x++)
    {
        for (int y=0; y<ysize; y++)
        {
            (*data)[x][y][0] = wall;
        }
    }

    for (int y=0; y<ysize; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            for (int z=0; z<zsize; z++)
            {
                (*data)[x][y][z].calcRho();
            }
        }
    }
}

//=========================== ACCESSORS ===========================

const DimSet3D Lattice3D::getSize()const
{
    DimSet3D pony = {{xsize, ysize,zsize}};
    return pony;
}

void Lattice3D::setData(const field3D& ndata, int x, int y, int z){
    data->resize(boost::extents[x][y][z]);
    *data = ndata;
    xsize = x;
    ysize = y;
    zsize = z;
}

void Lattice3D::setCell(int x, int y, int z, const Cell3D& ncell)
{
    if (y >= 0 && y < ysize && x >= 0 && x < xsize && z >= 0 && z < zsize) (*data)[x][y][z] = ncell;
}

void Lattice3D::setF(int x, int y, int z, int color, const array3D& nf)
{
    DistributionSetType3D f = (*data)[x][y][z].getF();
    f[color] = nf;
    (*data)[x][y][z].setF(f);
}

void Lattice3D::setF(int x, int y, int z, int color, int pos, double value)
{
    DistributionSetType3D f = (*data)[x][y][z].getF();
    f[color][pos] = value;
    (*data)[x][y][z].setF(f);
}

//=========================== OPERATOR ===========================

Lattice3D& Lattice3D::operator=(const Lattice3D& other){
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1], other.getSize()[2]);
    this->setParams(other.getParams());

    return *this;
}

const bool Lattice3D::operator==(const Lattice3D& other)const
{
    bool exit = true;
    DimSet3D extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1] && zsize == extent[2])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        field3D otherData = other.getData();
        for (int x = 0; x< xsize;x++)
        {
            for (int y = 0; y<ysize;y++)
            {
                for (int z = 0; z<zsize;z++)
                {
                    if (!((*data)[x][y][z]==otherData[x][y][z]))
                        {
                            exit = false;
                            break;
                        } // endif
                } // endfor z
            } // endfor y
        } //endfor x
    } // endif
    else exit = false;
   
    return exit;
}

// ///////////////////////////// PRIVATE /////////////////////////////

//=========================== OPERATIONS ===========================

void Lattice3D::linearIndex(int index, int& x, int& y, int& z)const
{
    x = (index)%xsize;
    y = ((index)/xsize)%ysize;
    z = (index)/(xsize * ysize);
}

void Lattice3D::streamAndBouncePull(Cell3D& tCell, const direction3D& dir)const
{
    const DistributionSetType3D f = tCell.getF();
    DistributionSetType3D ftmp;
    for (int color = 0; color<=1;color++)
    {
        ftmp[color][0] = (*data)[ dir[0].x ][ dir[0].y][ dir[0].z].getF()[color][0];

        for (int i=1;i<19;i++)
        {
            if ((*data)[ dir[PULL_INDEX_3D[i]].x][ dir[PULL_INDEX_3D[i]].y][ dir[PULL_INDEX_3D[i]].z].getIsSolid() == false)
            {
                ftmp[color][i] = (*data)[ dir[PULL_INDEX_3D[i]].x ][ dir[PULL_INDEX_3D[i]].y ][ dir[PULL_INDEX_3D[i]].z].getF()[color][i];    // if(neighbor not solid?) -> stream
            }
            else ftmp[color][i] = f[color][PULL_INDEX_3D[i]]; // else -> bounce back
        } // end for i
    } // end for color 
    tCell.setF(ftmp);
}

///////////////////////////// C-STYLE /////////////////////////////

//=========================== OPERATIONS ===========================

const DistributionSetType3D eqDistro(const ColSet& rho_k, const VeloSet3D& u, const DistributionSetType3D& phi)
{
    DistributionSetType3D feq;
    Vector3D velo = (u[0] * rho_k[0] + u[1] * rho_k[1]) * (1.0 / (rho_k[0]+rho_k[1])) ;
    double usqr = velo*velo;
   
    for (int i=0; i<19; i++)
    {
        double scal = velo*DIRECTION_3D[i];
        for (int color = 0; color<=1; color++)
        {
            feq[color][i] = rho_k[color] * ( phi[color][i] + WEIGHTS_3D[i] * ( 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
        }
    }
    return feq;
}

const DistributionSetType3D calculate_forcing_term(Vector3D G, VeloSet3D u)
{
    DistributionSetType3D forcing_term;
    for (int i=0; i<19; i++)
    {
        forcing_term[0][i] = 0;
        forcing_term[1][i] = WEIGHTS_3D[i] * (G * ( DIRECTION_3D[i] * (DIRECTION_3D[i] * u[1]) + DIRECTION_3D[i] - u[1])) ;
    }
    return forcing_term;
}