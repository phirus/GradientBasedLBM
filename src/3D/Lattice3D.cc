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

// //=========================== OPERATIONS ===========================

// void Lattice2D::equilibriumIni()
// {
//     Cell2D tmp;
//     DistributionSetType2D eqDis;

//     for (int j=0; j<ysize; j++)
//     {
//         for (int i=0; i<xsize; i++)
//         {
//             tmp = (*data)[i][j];
//             tmp.calcRho();
//             ColSet rho = tmp.getRho();
//             VeloSet2D u = tmp.getU();
//             eqDis = eqDistro(rho,u,param.getPhi());
//             tmp.setF(eqDis);
//             (*data)[i][j] = tmp;
//         }
//     }
// }

// void Lattice2D::balance(double& mass, double& momentum)const
// {
//     VeloSet2D u;
//     double rho;

//     mass = 0;
//     momentum = 0;

//     for (int j=0; j<ysize; j++)
//     {
//         for (int i=0; i<xsize; i++)
//         {
//             (*data)[i][j].calcRho();
//             rho = sum((*data)[i][j].getRho());
//             u = (*data)[i][j].getU();

//             mass += rho;
//             momentum += rho * sqrt(u[0]*u[0]);
//         }
//     }
// }

// void Lattice2D::mass_balance(double& liquid_mass, double& gas_mass)const
// {
//     ColSet rho;
//     liquid_mass = 0;
//     gas_mass = 0;

//     for (int j=0; j<ysize; j++)
//     {
//         for (int i=0; i<xsize; i++)
//         {
//             (*data)[i][j].calcRho();
//             rho = (*data)[i][j].getRho();
//             liquid_mass += rho[0];
//             gas_mass += rho[1];
//         }
//     }
// }

// void Lattice2D::overallRho()
// {
//     for (int j=0; j<ysize; j++)
//     {
//         for (int i=0; i<xsize; i++)
//         {
//             (*data)[i][j].calcRho();
//         }
//     }
// }

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

// const Vector2D Lattice2D::getGradient(int x, int y)const
// {
//     Vector2D grad(0,0);
//     double tmpDelta;

//     const direction2D dir = directions(x,y);
//     for (int q=0;q<13;q++)
//     {
//         tmpDelta = GRAD_WEIGHTS_2D[q] * (*data)[ dir[q].x ][ dir[q].y ].getDeltaRho();
//         grad.x += DIRECTION_2D[q].x * tmpDelta;
//         grad.y += DIRECTION_2D[q].y * tmpDelta;
//     }
//     return grad;
// }

// void Lattice2D::streamAll(int threads)
// {
//     field2D *newData = new field2D(boost::extents[xsize][ysize]);

//     omp_set_num_threads (threads);
//     const int range = xsize * ysize;

//     #pragma omp parallel
//     {
//         #pragma omp for        
//         for (int index = 0;  index < range; index++)
//         {
//             int x,y;
//             linearIndex(index,x,y);
//             const direction2D dir = directions(x,y);
//             Cell2D tmpCell = (*data)[x][y];

//             if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

//             tmpCell.calcRho();
//             #pragma omp critical(Zuweisung)
//             (*newData)[x][y] = tmpCell;
//         }
//     }
//     delete data;
//     data = newData;
// }

// bool Lattice2D::collideAll(int threads, bool gravity, bool isLimitActive)
// {
//     bool success(true);
//     field2D *newData = new field2D(boost::extents[xsize][ysize]);

//     omp_set_num_threads (threads);

//     const double beta = param.getBeta();
//     const DistributionSetType2D phi = param.getPhi();
//     const int range = xsize * ysize;
//     // const double rhoRedFixed = param.getRhoR();
//     const RelaxationPar2D relax = param.getRelaxation2D();
//     const double dt = param.getDeltaT();
//     // const double speedlimit = param.getSpeedlimit();

//     double g(0);
//     if(gravity == true) g = param.getG();
//     // if(gravity == true) g = 1e-1;


//     #pragma omp parallel
//     {
//         #pragma omp for
//         for (int index = 0;  index < range; index++)
//         {
//             int x,y;
//             linearIndex(index,x,y);
//             Cell2D tmpCell = (*data)[x][y];

//             if (tmpCell.getIsSolid() == false)
//             {
//                 DistributionSetType2D  fTmp;
//                 const DistributionSetType2D fCell = tmpCell.getF();

//                 const ColSet rho_k = tmpCell.getRho();
//                 const double rho = sum(rho_k);

//                 const Vector2D G(0 ,  g*(rho - rho_k[0]));
//                 // const Vector2D G(0 , - rho * g);

//                 VeloSet2D u = tmpCell.getU();

//                 if(gravity == true){
//                 u[0] = u[0] + G *  (dt/(2* rho)) ;
//                 u[1] = u[1] + G *  (dt/(2* rho)) ;
//                 }

//                 const DistributionSetType2D fEq = eqDistro(rho_k, u, phi);
//                 const DistributionSetType2D diff = distro_diff_2D(fCell, fEq);
            
//                 const double omega = param.getOmega(tmpCell.calcPsi());
//                 const Matrix2D relaxation_matrix(relax,omega);
                
//                 const Matrix2D forcing_factor = Matrix2D(true) - (relaxation_matrix*0.5);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
//                 const DistributionSetType2D first_forcing_term = forcing_factor * (TRAFO_MATRIX2D * calculate_forcing_term(G,u)); // F' = (I - 0.5 S) M F
//                 const DistributionSetType2D second_forcing_term = INV_TRAFO_MATRIX2D * first_forcing_term;    // M^{-1} F'

//                 const DistributionSetType2D single_phase_col = INV_TRAFO_MATRIX2D * (relaxation_matrix * (TRAFO_MATRIX2D * diff));
                          
//                 const ColSet A_k = param.getAk(omega);
//                 const Vector2D grad = getGradient(x,y);
//                 const double av = grad.Abs();

//                 double scal;
//                 double fges;
//                 double recolor;
//                 // double final_forcing_term(0);                

//                 for (int q=0; q<9; q++)
//                 {
//                     // gradient based two phase
//                     scal = grad*DIRECTION_2D[q];

//                     double gradient_collision(0);
//                     if (av > 0 ) gradient_collision = av/2 * (WEIGHTS_2D[q] * ( scal*scal )/(av*av) - B_2D[q]);

//                     for (int color=0;color<=1; color++)
//                     {
//                         // if (gravity == true) final_forcing_term = second_forcing_term[color][q] ;
//                         // fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q] + dt * final_forcing_term + A_k[color] * gradient_collision;
//                         fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];// + dt * final_forcing_term;
//                         if (gravity == true) fTmp[color][q] +=  dt * second_forcing_term[color][q];
//                     }
                    
//                     for (int color=0;color<=1; color++)
//                     {
//                         fTmp[color][q] += A_k[color] * gradient_collision;
//                     }

//                     fges = fTmp[0][q]+fTmp[1][q];

//                     // recoloring
//                     if (rho > 0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) *  grad.Angle(DIRECTION_2D[q])   * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q)); 
//                     else recolor = 0;

//                     if (rho > 0)
//                     {
//                         fTmp[0][q] = rho_k[0]/rho  * fges + recolor;
//                         fTmp[1][q] = rho_k[1]/rho * fges - recolor;
//                     }

//                     if (gravity == true){
//                         fTmp[0][q] +=  dt * second_forcing_term[0][q];
//                         fTmp[1][q] +=  dt * second_forcing_term[1][q];
//                     } 

//                     // if (fTmp[0][q] < 0) fTmp[0][q] = 0;
//                     // if (fTmp[1][q] < 0) fTmp[1][q] = 0;

//                 } // end for 

//                 tmpCell.setF(fTmp);
//                 tmpCell.calcRho();
//             }
//             #pragma omp critical(Zuweisung2)
//             (*newData)[x][y] = tmpCell;
//         }
//     }
//     if(success == true){
//         delete data;
//         data = newData;
//     }
//     else{
//         delete newData;
//     }     

//     return success;
// }

// void Lattice2D::closedBox()
// {
//     const Cell2D wall(0,0,true);

//     for (int x=0; x<xsize; x++)
//     {
//         (*data)[x][0] = wall;
//         (*data)[x][ysize-1] = wall;
//     }
//     for (int y=0; y<ysize; y++)
//     {
//         (*data)[0][y] = wall;
//         (*data)[xsize-1][y] = wall;
//     }

//     for (int y=0; y<ysize; y++)
//     {
//         for (int x=0; x<xsize; x++)
//         {
//             (*data)[x][y].calcRho();
//         }
//     }
// }

// void Lattice2D::bottomWall()
// {
//     const Cell2D wall(0,0,true);

//     for (int x=0; x<xsize; x++)
//     {
//         (*data)[x][0] = wall;
//     }

//     for (int y=0; y<ysize; y++)
//     {
//         for (int x=0; x<xsize; x++)
//         {
//             (*data)[x][y].calcRho();
//         }
//     }
// }

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

// void Lattice2D::streamAndBouncePull(Cell2D& tCell, const direction2D& dir)const
// {
//     const DistributionSetType2D f = tCell.getF();
//     DistributionSetType2D ftmp;
//     for (int color = 0; color<=1;color++)
//     {
//         ftmp[color][0] = (*data)[ dir[0].x ][ dir[0].y ].getF()[color][0];

//         if ((*data)[ dir[5].x ][ dir[5].y ].getIsSolid() == false)
//         {
//             ftmp[color][1] = (*data)[ dir[5].x ][ dir[5].y ].getF()[color][1];    // if(neighbor not solid?) -> stream
//         }
//         else ftmp[color][1] = f[color][5];                                                                                            // else -> bounce back

//         if ((*data)[ dir[6].x ][ dir[6].y ].getIsSolid() == false)
//         {
//             ftmp[color][2] = (*data)[ dir[6].x ][ dir[6].y ].getF()[color][2];
//         }
//         else ftmp[color][2] = f[color][6];

//         if ((*data)[ dir[7].x ][ dir[7].y ].getIsSolid() == false)
//         {
//             ftmp[color][3] = (*data)[ dir[7].x ][ dir[7].y ].getF()[color][3];
//         }
//         else ftmp[color][3] = f[color][7];

//         if ((*data)[ dir[8].x ][ dir[8].y ].getIsSolid() == false)
//         {
//             ftmp[color][4] = (*data)[ dir[8].x ][ dir[8].y ].getF()[color][4];
//         }
//         else ftmp[color][4] = f[color][8];

//         if ((*data)[ dir[1].x ][ dir[1].y ].getIsSolid() == false)
//         {
//             ftmp[color][5] = (*data)[ dir[1].x ][ dir[1].y ].getF()[color][5];
//         }
//         else ftmp[color][5] = f[color][1];

//         if ((*data)[ dir[2].x ][ dir[2].y ].getIsSolid() == false)
//         {
//             ftmp[color][6] = (*data)[ dir[2].x ][ dir[2].y ].getF()[color][6];
//         }
//         else ftmp[color][6] = f[color][2];

//         if ((*data)[ dir[3].x ][ dir[3].y ].getIsSolid() == false)
//         {
//             ftmp[color][7] = (*data)[ dir[3].x ][ dir[3].y ].getF()[color][7];
//         }
//         else ftmp[color][7] = f[color][3];

//         if ((*data)[ dir[4].x ][ dir[4].y ].getIsSolid() == false)
//         {
//             ftmp[color][8] = (*data)[ dir[4].x ][ dir[4].y ].getF()[color][8];
//         }
//         else ftmp[color][8] = f[color][4];
//     }
//     tCell.setF(ftmp);
// }

// ///////////////////////////// C-STYLE /////////////////////////////

// //=========================== OPERATIONS ===========================

// const DistributionSetType2D eqDistro(const ColSet& rho_k, const VeloSet2D& u, const DistributionSetType2D& phi)
// {
//     DistributionSetType2D feq;
//     Vector2D velo = (u[0] * rho_k[0] + u[1] * rho_k[1]) * (1.0 / (rho_k[0]+rho_k[1])) ;
//     double usqr = velo*velo;
   
//     for (int i=0; i<9; i++)
//     {
//         double scal = velo*DIRECTION_2D[i];
//         for (int color = 0; color<=1; color++)
//         {
//             feq[color][i] = rho_k[color] * ( phi[color][i] + WEIGHTS_2D[i] * ( 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr));
//         }
//     }
//     return feq;
// }

// const DistributionSetType2D calculate_forcing_term(Vector2D G, VeloSet2D u)
// {
//     DistributionSetType2D forcing_term;
//     for (int i=0; i<9; i++)
//     {
//         forcing_term[0][i] = 0;
//         forcing_term[1][i] = WEIGHTS_2D[i] * (G * ( DIRECTION_2D[i] * (DIRECTION_2D[i] * u[1]) + DIRECTION_2D[i] - u[1])) ;
//     }
//     return forcing_term;
// }