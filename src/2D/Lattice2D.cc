#include"Lattice2D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice2D::Lattice2D(int x_size, int y_size,double fzero_dense, double fzero_dilute):
xsize(x_size)
,ysize(y_size)
,data(new field2D(boost::extents[xsize][ysize]))
,param()
,bound()
,bubblebox(0,0,xsize,ysize)
,offset(0)
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
,bound(other.getBoundaries())
,bubblebox(other.getBubbleBox())
,offset(other.getOffset())
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
//        tmpDelta = GRAD_WEIGHTS_2D[q] * (*data)[ dir[q].x ][ dir[q].y ].getDeltaRho();
        tmpDelta = GRAD_WEIGHTS_2D[q] * (*data)[ dir[q].x ][ dir[q].y ].calcPsi();
        grad.x += DIRECTION_2D[q].x * tmpDelta;
        grad.y += DIRECTION_2D[q].y * tmpDelta;
    }
    return grad;
}

const Vector2D Lattice2D::getSurfaceNormal(int x, int y) const
{
    const Vector2D grad = getGradient(x,y);
    const double abs = grad.Abs();
    Vector2D n(0,0);
    if (abs > 1e-3) n = grad * (1.0 / abs);
    return n;
}

void Lattice2D::streamAll(int threads)
{
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);
    //const int range = xsize * ysize;

    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(dynamic, 100)
        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                const direction2D dir = directions(x,y);
                Cell2D tmpCell = (*data)[x][y];

                if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);

                tmpCell.calcRho();
                (*newData)[x][y] = tmpCell;
            }
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
    const double dt = param.getDeltaT();

    double g(0);
    if(gravity == true) g = param.getG();

    const Matrix2D relaxation_matrix_outside(param.getRelaxation2D(1),false);
    const Matrix2D single_relax_outside = INV_TRAFO_MATRIX2D * relaxation_matrix_outside * TRAFO_MATRIX2D;

    #pragma omp parallel firstprivate(single_relax_outside)
    {
        #pragma omp for collapse(2) schedule(dynamic, 100)
        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                Cell2D tmpCell = (*data)[x][y];

                if (tmpCell.getIsSolid() == false && isBoundary(x, y) == false)
                {
                    //const direction2D dir = directions(x,y);

                    DistributionSetType2D  fTmp;
                    const DistributionSetType2D fCell = tmpCell.getF();

                    const ColSet rho_k = tmpCell.getRho();
                    const double rho = sum(rho_k);
                    VeloSet2D u = tmpCell.getU();

                    double fges;
 
                    if(bubblebox.isInBox(x,y) == false)
                    {
                        const DistributionSetType2D fEq = eqDistro(rho_k, u, phi);
                        const DistributionSetType2D diff = distro_diff_2D(fCell, fEq);
                        const DistributionSetType2D single_phase_col =  single_relax_outside * diff;                        

                        for (int q=0; q<9; q++)
                        {
                            for (int color=0;color<=1; color++)
                            {

                                fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];
                            }
                            fges = fTmp[0][q]+fTmp[1][q];

                            if (rho > 0)
                            {
                                fTmp[0][q] = rho_k[0]/rho  * fges;
                                fTmp[1][q] = rho_k[1]/rho * fges;
                            }
                        } // end for
                    } // end if
                    else
                    {
                        const double psi = tmpCell.calcPsi();
                        const RelaxationPar2D relax = param.getRelaxation2D(psi);
                        const double omega = param.getOmega(psi);
                        
                        const Matrix2D relaxation_matrix(relax,false);
                        const Matrix2D single_relax = INV_TRAFO_MATRIX2D * relaxation_matrix * TRAFO_MATRIX2D;
                        const Matrix2D forcing_factor(relax,true);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
                        const Matrix2D forcing_relax = INV_TRAFO_MATRIX2D * forcing_factor * TRAFO_MATRIX2D; 
                        
                        const ColSet A_k = param.getAk(omega);

                        const Vector2D G(0 ,  g*(rho - rho_k[0]));
                        
                        if(gravity == true)
                        {
                            u[0] = u[0] + G *  (dt/(2* rho)) ;
                            u[1] = u[1] + G *  (dt/(2* rho)) ;
                        }

                        const DistributionSetType2D fEq = eqDistro(rho_k, u, phi);
                        const DistributionSetType2D diff = distro_diff_2D(fCell, fEq);
                        const DistributionSetType2D single_phase_col =  single_relax * diff;

                        const DistributionSetType2D second_forcing_term = forcing_relax * calculate_forcing_term(G,u);    // M^{-1} F'
                        const Vector2D grad = getGradient(x,y);
                        const double av = grad.Abs();

                        double scal;
                        double recolor;

                        // q = 0
                        for (int color=0;color<=1; color++)
                        {
                            fTmp[color][0] =  fCell[color][0] - single_phase_col[color][0]; //single
                            if (gravity == true) fTmp[color][0] +=  dt * second_forcing_term[color][0]; // forcing term
                            fTmp[color][0] += A_k[color] * (- (av/2 * B_2D[0])); // perturbation
                        }
                        fges = fTmp[0][0]+fTmp[1][0];
                        if (rho > 0)
                        {
                            fTmp[0][0] = rho_k[0]/rho  * fges;
                            fTmp[1][0] = rho_k[1]/rho * fges;
                        }

                        for (int q=1; q<9; q++)
                        {
                            scal = grad*DIRECTION_2D[q];
                            double gradient_collision(0);
                            if (av > 0 ) gradient_collision = av/2 * (WEIGHTS_2D[q] * ( scal*scal )/(av*av) - B_2D[q]);

                            for (int color=0;color<=1; color++)
                            {
                                fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];
                                if (gravity == true) fTmp[color][q] +=  dt * second_forcing_term[color][q];
                                fTmp[color][q] += A_k[color] * gradient_collision;
                            }

                            fges = fTmp[0][q]+fTmp[1][q];

                            // recoloring
                            if (rho > 0 && av>0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) * (scal /(av * DIRECTION_ABS_2D[q])) * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q)); 
                            else recolor = 0;

                            if (rho > 0)
                            {
                                fTmp[0][q] = rho_k[0]/rho  * fges + recolor;
                                fTmp[1][q] = rho_k[1]/rho * fges - recolor;
                            }

                        } // end for
                    } // end else

                    tmpCell.setF(fTmp);
                }
                (*newData)[x][y] = tmpCell;
            }
        }
    }
    
    
    if(success == true)
    {
        delete data;
        data = newData;
    }
    else
    {
        delete newData;
    }     

    return success;
}

void Lattice2D::evaluateBoundaries(int threads)
{
    omp_set_num_threads (threads);

    #pragma omp parallel
    {
        // north boundary
        if(bound.north.getType() > 1)
        {
            int lowerX = 0;
            int upperX = xsize;

            if(bound.west.getType()  > 1) // north west corner
            {
                (*data)[0][ysize-1] = cornerNorthWest((*data)[0][ysize-1], bound.north.getRho());
                lowerX = 1;
            }
            else if(bound.west.getType() == bounceback)
            {
                (*data)[1][ysize-1] = cornerNorthWest((*data)[1][ysize-1], bound.north.getRho());
                lowerX = 2;
            }

            if(bound.east.getType() > 1)
            {
                (*data)[xsize-1][ysize-1] = cornerNorthEast((*data)[xsize-1][ysize-1], bound.north.getRho());
                upperX = xsize - 1;
            }
            else if(bound.east.getType() == bounceback)
            {
                (*data)[xsize-2][ysize-1] = cornerNorthEast((*data)[xsize-2][ysize-1], bound.north.getRho());
                upperX = xsize - 2;
            }

            if(bound.north.getType() == pressure)
            {
                ColSet rho = bound.north.getRho();
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    (*data)[x][ysize-1] = boundaryNorthPres((*data)[x][ysize-1], rho);
                }
            }
            if(bound.north.getType() == velocity)
            {
                double u_y = bound.north.getVelocity()[0].y;
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    (*data)[x][ysize-1] = boundaryNorthVelo((*data)[x][ysize-1], u_y);
                }
            }
            if(bound.north.getType() == shear)
            {
                double u_0 = bound.west.getVelocity()[0].y;
                double u_1 = bound.east.getVelocity()[0].y;

                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    double u_y = get_shearvelocity_x(x,u_0,u_1);
                    (*data)[x][ysize-1] = boundaryNorthVelo((*data)[x][ysize-1], u_y);
                }
            }
        } // end north

        // south boundary
        if(bound.south.getType() > 1)
        {
            int lowerX = 0;
            int upperX = xsize;
    
            if(bound.west.getType() > 1)
            {
                (*data)[0][0] = cornerSouthWest((*data)[0][0], bound.south.getRho());
                lowerX = 1;
            }
            else if (bound.west.getType() == bounceback)
            {
                (*data)[1][0] = cornerSouthWest((*data)[1][0], bound.south.getRho());
                lowerX = 2;
            }
    
            if(bound.east.getType() > 1)
            {
                (*data)[xsize-1][0] = cornerSouthEast((*data)[xsize-1][0], bound.south.getRho());
                upperX = xsize - 1;
            }
            else if(bound.east.getType() == bounceback)
            {
                (*data)[xsize-2][0] = cornerSouthEast((*data)[xsize-2][0], bound.south.getRho());
                upperX = xsize - 2;
            }

            if(bound.south.getType() == pressure)
            {
                ColSet rho = bound.south.getRho();
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    (*data)[x][0] = boundarySouthPres((*data)[x][0],rho);
                }
            }

            if(bound.south.getType() == velocity)
            {
                double u_y = bound.south.getVelocity()[0].y;
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    (*data)[x][0] = boundarySouthVelo((*data)[x][0],u_y);
                }
            }
            if(bound.south.getType() == shear)
            {
                double u_0 = bound.west.getVelocity()[0].y;
                double u_1 = bound.east.getVelocity()[0].y;
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    double u_y = get_shearvelocity_x(x,u_0,u_1);
                    (*data)[x][0] = boundarySouthVelo((*data)[x][0],u_y);
                }
            }
        } // end south 

        // west boundary
        if(bound.west.getType() > 1)
        {
            int lowerY = 0;
            int upperY = ysize;

            if (bound.south.getType() > 1) lowerY = 1;
            else if (bound.south.getType() == bounceback)
            {
                (*data)[0][1] = cornerSouthWest((*data)[0][1], bound.west.getRho());
                lowerY = 2;
            }

            if (bound.north.getType() > 1) upperY = ysize - 1;
            else if (bound.north.getType() == bounceback)
            {
                (*data)[0][ysize-2] = cornerNorthWest((*data)[0][ysize-2], bound.west.getRho());
                upperY = ysize - 2;
            }

            if(bound.west.getType() == pressure)
            {
                ColSet rho = bound.west.getRho();
                #pragma omp for schedule(static,10) nowait
                for (int y=lowerY; y< upperY; y++)
                {
                    (*data)[0][y] = boundaryWestPres((*data)[0][y],rho);
                }
            }
            if(bound.west.getType() == velocity)
            {
                double u_x = bound.west.getVelocity()[0].x;
                #pragma omp for schedule(static,10) nowait
                for (int y=lowerY; y< upperY; y++)
                {
                    (*data)[0][y] = boundaryWestVelo((*data)[0][y], u_x);
                }
            }    
        } // end west
    
        // east boundary
        if(bound.east.getType() > 1)
        {
            int lowerY = 0;
            int upperY = ysize;
    
            if (bound.south.getType() > 1) lowerY = 1;
            else if (bound.south.getType() == bounceback)
            {
                (*data)[xsize-1][1] = cornerSouthEast((*data)[xsize-1][1], bound.east.getRho());
                lowerY = 2;
            }

            if (bound.north.getType() > 1) upperY = ysize - 1;
            else if (bound.north.getType() == bounceback)
            {
                (*data)[xsize-1][ysize-2] = cornerNorthEast((*data)[xsize-1][ysize-2], bound.east.getRho());
                upperY = ysize - 2;
            }

            if(bound.east.getType() == pressure)
            {
                #pragma omp for schedule(static,10) nowait
                for (int y=lowerY; y< upperY; y++)
                {
                    ColSet rho = bound.east.getRho();
                    (*data)[xsize-1][y] = boundaryEastPres((*data)[xsize-1][y], rho);
                }
            }
            if(bound.east.getType() == velocity)
            {
                double u_x = bound.east.getVelocity()[0].x;
                #pragma omp for schedule(static,10) nowait
               for (int y=lowerY; y< upperY; y++)
                {
                    (*data)[xsize-1][y] = boundaryEastVelo((*data)[xsize-1][y], u_x) ;
                }
            }           
        } //end east
    }
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

void Lattice2D::setShearProfile(double gradient, double offset)
{
    const int range = xsize * ysize;
    const double m = gradient;
    const double n = offset; // - (m * ( xsize *  param.getDeltaX() ) );
    for (int index = 0;  index < range; index++)
        {
            int x,y;
            linearIndex(index,x,y);
            Cell2D tmpCell = (*data)[x][y];
            tmpCell.calcRho();
            VeloSet2D u;
            u[0] = Vector2D(0,m*x + n);
            u[1] = Vector2D(0,0);
            //const Vector2D u(0 , m*x + n);   
            const ColSet rho = tmpCell.getRho();
            const DistributionSetType2D phi = param.getPhi2D();

            tmpCell.setF(eqDistro(rho,u, phi));
            (*data)[x][y] = tmpCell;
        }

}

//=========================== ACCESSORS ===========================

const DimSet2D Lattice2D::getSize()const
{
    DimSet2D pony = {{xsize, ysize}}; 
    return pony;
}

const field2D Lattice2D::getData(int cutoff)const
{
    field2D newData(boost::extents[xsize][ysize - cutoff]);

    for (int y=0; y<ysize - cutoff; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            newData[x][y] = (*data)[x][y + cutoff];
        }
    }

    return newData;
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

void Lattice2D::setBoundaries(const Boundaries& newBound)
{
    bound = newBound;
    buildWalls();
}

void Lattice2D::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}

//=========================== LATTICE CUTOUT ===========================

const Lattice2D Lattice2D::latticeCutOff(int cutoff)const
{
    const field2D nd = getData(cutoff);
    BubbleBox2D newBB = bubblebox;
    newBB.setY(bubblebox.getY() - cutoff);

    Lattice2D newLattice;
    newLattice.setData(nd,xsize,ysize - cutoff);
    newLattice.setParams(param);
    newLattice.setBoundaries(bound);
    newLattice.setBubbleBox(newBB);
    newLattice.setOffset(offset + cutoff);

    return newLattice;
}

const Lattice2D Lattice2D::latticeAppend(int append, double rho_red, double rho_blue)const
{

    field2D nd(boost::extents[xsize][ysize + append]);

    Cell2D newCell(eqDistro({{rho_red,rho_blue}}, VeloSet2D(), param.getPhi2D()));
    newCell.calcRho();

    for (int y=0; y<ysize + append; y++)
    {
        for (int x=0; x<xsize; x++)
        {
            if(y < ysize) nd[x][y] = (*data)[x][y];
            else nd[x][y] = newCell;;
        }
    }

    Lattice2D newLattice;
    newLattice.setData(nd,xsize,ysize + append);
    newLattice.setParams(param);
    newLattice.setBoundaries(bound);
    newLattice.setBubbleBox(bubblebox);
    newLattice.setOffset(offset);

    return newLattice;
}


const std::vector<int> Lattice2D::findBubbleCells()const
{
    std::vector<int> indices;
    const int range = xsize * ysize;
    int x,y;

    for (int index = 0;  index < range; index++)
    {
        linearIndex(index,x,y);

        const Cell2D tmpCell = (*data)[x][y];

        if (tmpCell.calcPsi() < BUBBLE_CRITERION && tmpCell.getIsSolid() == false)
        {
            indices.push_back(index);
        }        
    }
    return indices;
}

void Lattice2D::copyCellsFromOther(const Lattice2D& other, const std::vector<int>& indices)
{
    int x,y;
    for (int index : indices)
    {
        linearIndex(index,x,y);
        (*data)[x][y] = other.getCell(x,y);
    }
}

const boost::array<Vector2D,3> Lattice2D::getBubbleData()const
{
    Vector2D momentum(0,0);
    Vector2D tmp_position(0,0);
    Vector2D force(0,0);
    double rho_sum(0);


    for (int x = 0; x<xsize;x++)
    {
        for (int y = 0; y<ysize;y++)
        {
            const Cell2D tmp_cell = (*data)[x][y];
            const VeloSet2D tmp_velo = tmp_cell.getU();
            const ColSet tmp_rho = tmp_cell.getRho();
            if ( tmp_cell.calcPsi() < 0.99 && tmp_cell.getIsSolid() == false) 
            {
                momentum = momentum + (tmp_velo[1] * tmp_rho[1]);
                tmp_position = tmp_position + (Vector2D(x,y) * tmp_rho[1]);
                force = force +  getSurfaceNormal(x, y) * getDivergence(x, y);
                rho_sum += tmp_rho[1];
            }
        }
    }

    const Vector2D velocity = momentum * (1.0/ rho_sum);
    const Vector2D position = tmp_position * (1.0/ rho_sum);
    
    boost::array<Vector2D,3> result = {{position,velocity,force}};
    return result;

}

const Vector2D Lattice2D::getPGesN(int x, int y)const
{
    const Cell2D tmp = (*data)[x][y];
    const double p_ges = tmp.getPressureTensor(0).getTrace() + tmp.getPressureTensor(1).getTrace();
    const Vector2D n = getSurfaceNormal(x, y) * (-1);

    return n * p_ges;
}

const Vector2D Lattice2D::getP1N(int x, int y)const
{
    const Cell2D tmp = (*data)[x][y];
    double p_ges = tmp.getPressureTensor(0).getTrace();
    Vector2D n = getSurfaceNormal(x, y) * (-1);

    return n * p_ges;
}


const double Lattice2D::getDivergence(int x, int y)const
{
    double div = 0;
    double p_x,p_y;

    const direction2D dir = directions(x,y);
    for (int q=0;q<13;q++)
    {
        Vector2D tmpP = getP1N(dir[q].x, dir[q].y);
        p_x = GRAD_WEIGHTS_2D[q] * tmpP.x;
        p_y = GRAD_WEIGHTS_2D[q] * tmpP.y;
        div += DIRECTION_2D[q].x * p_x + DIRECTION_2D[q].y * p_y;
    }
    return div;
}

const Vector2D Lattice2D::getResultingForce()const
{
    Vector2D force(0,0);
    const std::vector<int> indices = findBubbleCells();

    int x,y;

    for(int index : indices)
    {
        linearIndex(index, x, y);
        force = force +  getSurfaceNormal(x, y) * getDivergence(x, y) ;
    }
    return force;
}

//====================== BOUNDARY TREATMENT ======================

const bool Lattice2D::isBoundary(int x, int y)const
{
    return  (x == 0 && bound.west.getType()!= periodic) || (x == xsize-1 && bound.east.getType()!= periodic) || (y == ysize-1 && bound.north.getType()!= periodic) || (y == 0 && bound.south.getType()!= periodic);
}

void Lattice2D::buildWalls()
{
    Cell2D wall(0,0,true);

    if(bound.north.getType() == bounceback)
    {
        Vector3D u3D = bound.north.getVelocity()[0];
        Vector2D u_wall(u3D.x,u3D.y);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            (*data)[x][ysize-1] = wall;
        }
    }

    if(bound.south.getType() == bounceback)
    {
        Vector3D u3D = bound.south.getVelocity()[0];
        Vector2D u_wall(u3D.x,u3D.y);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            (*data)[x][0] = wall;
        }
    }

    if(bound.west.getType() == bounceback)
    {
        Vector3D u3D = bound.west.getVelocity()[0];
        Vector2D u_wall(u3D.x,u3D.y);
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            (*data)[0][y] = wall;
        }
    }

    if(bound.east.getType() == bounceback)
    {
        Vector3D u3D = bound.east.getVelocity()[0];
        Vector2D u_wall(u3D.x,u3D.y);
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            (*data)[xsize-1][y] = wall;
        }
    }    

}

const Cell2D Lattice2D::boundaryNorthPres(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    ColSet u_y;
    
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0)
        {
            u_y[color] = -1.0 + (f[color][0] + f[color][1] + f[color][5] + 2* (f[color][2] + f[color][3] + f[color][4])) / rho[color];
            f[color][7] = f[color][3] - 2.0/3.0 * rho[color]*u_y[color];
            f[color][6] = - rho[color]*u_y[color]/6.0 + f[color][2] + (f[color][1] - f[color][5])/2.0;
            f[color][8] = - rho[color]*u_y[color]/6.0 + f[color][4] + (f[color][5] - f[color][1])/2.0;
        }
        else
        {
            u_y[color] = 0;
            f[color][7] = 0;
            f[color][6] = 0;
            f[color][8] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::boundaryNorthVelo(const Cell2D& tmp, double uy)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(0,uy),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D::boundarySouthPres(const Cell2D& tmp, ColSet rho)const
{        
    DistributionSetType2D f = tmp.getF();
    ColSet u_y;
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0 )
        {
            u_y[color] = 1.0 - (f[color][0] + f[color][1] + f[color][5] + 2* (f[color][6] + f[color][7] + f[color][8])) / rho[color];
            f[color][3] = f[color][7] + 2.0/3.0 * rho[color]*u_y[color];
            f[color][2] = rho[color]*u_y[color]/6.0 + f[color][6] + (f[color][5] - f[color][1])/2.0;
            f[color][4] = rho[color]*u_y[color]/6.0 + f[color][8] + (f[color][1] - f[color][5])/2.0;
        }
        else
        {
            u_y[color] = 0;
            f[color][3] = 0;
            f[color][2] = 0;
            f[color][4] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::boundarySouthVelo(const Cell2D& tmp, double uy)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(0,uy),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D::boundaryWestPres(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    ColSet u_x;    
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0)
        {
            u_x[color] = 1.0 - (f[color][0] + f[color][3] + f[color][7] + 2* (f[color][4] + f[color][5] + f[color][6])) / rho[color];
            f[color][1] = f[color][5] + 2.0/3.0 * rho[color]*u_x[color];
            f[color][2] = rho[color]*u_x[color]/6.0 + f[color][6] + (f[color][7] - f[color][3])/2.0;
            f[color][8] = rho[color]*u_x[color]/6.0 + f[color][4] + (f[color][3] - f[color][7])/2.0;
        }
        else
        {
            u_x[color] = 0;
            f[color][1] = 0;
            f[color][2] = 0;
            f[color][8] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::boundaryWestVelo(const Cell2D& tmp, double ux)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(ux,0),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D::boundaryEastPres(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    ColSet u_x;
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0)
        {
            u_x[color] = -1.0 + (f[color][0] + f[color][3] + f[color][7] + 2* (f[color][1] + f[color][2] + f[color][8])) / rho[color];
            f[color][5] = f[color][1] - 2.0/3.0 * rho[color]*u_x[color];
            f[color][6] = - rho[color]*u_x[color]/6.0 + f[color][2] + (f[color][3] - f[color][7])/2.0;
            f[color][4] = - rho[color]*u_x[color]/6.0 + f[color][8] + (f[color][7] - f[color][3])/2.0;
        }
        else
        {
            u_x[color] = 0;
            f[color][5] = 0;
            f[color][6] = 0;
            f[color][4] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::boundaryEastVelo(const Cell2D& tmp, double ux)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(ux,0),Vector2D(0,0)}}, param.getPhi2D()));
}

// corners
const Cell2D Lattice2D::cornerNorthWest(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();

    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0 )
        {
            f[color][1] = f[color][5];
            f[color][7] = f[color][3];
            f[color][8] = f[color][4];
            f[color][2] = 0.5 * (rho[color] - f[color][0] )  - f[color][3] - f[color][4] - f[color][5];
            f[color][6] = f[color][2];
        }
        else
        {
            f[color][1] = 0;
            f[color][2] = 0;
            f[color][6] = 0;
            f[color][7] = 0;
            f[color][8] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::cornerNorthEast(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0 )
        {
            f[color][5] = f[color][1];
            f[color][6] = f[color][2];
            f[color][7] = f[color][3];
            f[color][4] = 0.5 * (rho[color] - f[color][0] )  - f[color][1] - f[color][2] - f[color][3];
            f[color][8] = f[color][4];
        }
        else
        {
            f[color][4] = 0;
            f[color][5] = 0;
            f[color][6] = 0;
            f[color][7] = 0;
            f[color][8] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::cornerSouthWest(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0 )
        {
            f[color][1] = f[color][5];
            f[color][2] = f[color][6];
            f[color][3] = f[color][7];
            f[color][4] = 0.5 * (rho[color] - f[color][0] )  - f[color][5] - f[color][6] - f[color][7];
            f[color][8] = f[color][4];
        }
        else
        {
            f[color][1] = 0;
            f[color][2] = 0;
            f[color][3] = 0;
            f[color][4] = 0;
            f[color][8] = 0;
        }
    }
    return Cell2D(f);
}

const Cell2D Lattice2D::cornerSouthEast(const Cell2D& tmp, ColSet rho)const
{
    DistributionSetType2D f = tmp.getF();
    for(int color = 0; color <= 1; color++)
    {
        if(rho[color]>0 )
        {
            f[color][3] = f[color][7];
            f[color][4] = f[color][8];
            f[color][5] = f[color][1];
            f[color][2] = 0.5 * (rho[color] - f[color][0] )  - f[color][1] - f[color][7] - f[color][8];
            f[color][6] = f[color][2];
        }
        else
        {
            f[color][3] = 0;
            f[color][4] = 0;
            f[color][5] = 0;
            f[color][2] = 0;
            f[color][6] = 0;
        }
    }
    return Cell2D(f);
}


//=========================== OPERATOR ===========================

Lattice2D& Lattice2D::operator=(const Lattice2D& other)
{
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1]);
    this->setParams(other.getParams());
    this->setBoundaries(other.getBoundaries());
    this->setBubbleBox(other.getBubbleBox());
    this->setOffset(other.getOffset());

    return *this;
}

const bool Lattice2D::operator==(const Lattice2D& other)const
{
    bool exit = true;
    DimSet2D extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        Boundaries bOther = other.getBoundaries();
        if (!(bound == bOther)) exit = false;

        if (!(bubblebox == other.getBubbleBox())) exit = false;

        if (!(offset == other.getOffset())) exit = false;

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
                ftmp[color][i] = f[color][PULL_INDEX_2D[i]] + (2.0 * 3.0 * WEIGHTS_2D[i] * rho[color] * (DIRECTION_2D[i] * neighbor.getU()[0]) ) ;
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
        forcing_term[1][i] = WEIGHTS_2D[i] * (G * ( DIRECTION_2D[i] * (DIRECTION_2D[i] * u[1]) * 9 + (DIRECTION_2D[i] - u[1]) * 3)) ;
    }
    return forcing_term;
}

// const Vector2D getGradientFun(const direction2D& dir, const field2D& data)
// {
//     Vector2D grad(0,0);
//     double tmpDelta;

//     for (int q=0;q<13;q++)
//     {
//         tmpDelta = GRAD_WEIGHTS_2D[q] * data[ dir[q].x ][ dir[q].y ].getDeltaRho();
//         grad = grad + (DIRECTION_2D[q] * tmpDelta);
//     }
//     return grad;
// }
