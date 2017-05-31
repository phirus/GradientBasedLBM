#include"Lattice2D_no_inter.h"

// ///////////////////////////// PUBLIC /////////////////////////////

// //=========================== LIFECYCLE ===========================

Lattice2D_no_inter::Lattice2D_no_inter(int x_size, int y_size,double fzero_dense, double fzero_dilute):
xsize(x_size)
,ysize(y_size)
,data(new field2D(boost::extents[xsize][ysize]))
,param()
,bound()
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

Lattice2D_no_inter::Lattice2D_no_inter(const Lattice2D_no_inter& other):
xsize(other.getSize()[0])
,ysize(other.getSize()[1])
,data(new field2D(boost::extents[xsize][ysize]))
,param(other.getParams())
,bound(other.getBoundaries())
,offset(other.getOffset())
{
    (*data) = other.getData();
}


Lattice2D_no_inter::~Lattice2D_no_inter(){
    delete data;
    data = NULL;
}

// //=========================== OPERATIONS ===========================

void Lattice2D_no_inter::balance(double& mass, double& momentum)const
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

void Lattice2D_no_inter::mass_balance(double& liquid_mass, double& gas_mass)const
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

void Lattice2D_no_inter::overallRho()
{
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            (*data)[i][j].calcRho();
        }
    }
}

direction2D Lattice2D_no_inter::directions(int x, int y)const
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

void Lattice2D_no_inter::streamAll(int threads)
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

bool Lattice2D_no_inter::collideAll(int threads)
{
    bool success(true);
    field2D *newData = new field2D(boost::extents[xsize][ysize]);

    omp_set_num_threads (threads);

//    const Matrix2D relaxation_matrix(param.getRelaxation2D(1),false);
//    const Matrix2D single_relax = INV_TRAFO_MATRIX2D * relaxation_matrix * TRAFO_MATRIX2D;

    #pragma omp parallel // firstprivate(single_relax)
    {
        #pragma omp for collapse(2) schedule(dynamic, 100)
        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                Cell2D tmpCell = (*data)[x][y];

                if (tmpCell.getIsSolid() == false && isBoundary(x, y) == false)
                {
                    //const DistributionSetType2D fCell = tmpCell.getF();

                    const ColSet rho_k = tmpCell.getRho();
                    const double rho = sum(rho_k);
                    VeloSet2D u_k = tmpCell.getU();

                    DistributionSetType2D fNew;
                    if(rho>0)
                    {
                        Vector2D u = (u_k[0] * rho_k[0] + u_k[1] * rho_k[1]) * (1.0/rho);
                        const array2D fEq = eqDistro_no_inter(rho, u);
                        // const array2D fGes = array_add_2D(fCell[0],fCell[1]);
                        // const array2D diff = array_diff_2D(fGes, fEq);
                        // const array2D single_phase_col =  single_relax * diff;
                        // const array2D  fTmp = array_diff_2D(fTmp, single_phase_col);
                        fNew[0] = array_times_2D(fEq, rho_k[0]/rho);
                        fNew[1] = array_times_2D(fEq, rho_k[1]/rho);
                    }
                    tmpCell.setF(fNew);
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

void Lattice2D_no_inter::evaluateBoundaries(int threads)
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

void Lattice2D_no_inter::closedBox()
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

void Lattice2D_no_inter::bottomWall()
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

void Lattice2D_no_inter::genericWall(std::vector<double> x, std::vector<double> y, const Vector2D& u_w)
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

void Lattice2D_no_inter::lidDrivenCavity(const Vector2D& u_w)
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

void Lattice2D_no_inter::shearWall(const Vector2D& u_w)
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

void Lattice2D_no_inter::setShearProfile(double gradient, double offset)
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

const DimSet2D Lattice2D_no_inter::getSize()const
{
    DimSet2D pony = {{xsize, ysize}}; 
    return pony;
}

const field2D Lattice2D_no_inter::getData(int cutoff)const
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

const Lattice2D Lattice2D_no_inter::getLattice2D()const
{
    Lattice2D l(xsize,ysize);
    l.setData(getData(),xsize, ysize);
    l.setParams(getParams());
    l.setBoundaries(getBoundaries());

    return l;
}

void Lattice2D_no_inter::setData(const field2D& ndata, int x, int y){
    data->resize(boost::extents[x][y]);
    *data = ndata;
    xsize = x;
    ysize = y;
}

void Lattice2D_no_inter::setCell(int x, int y, const Cell2D& ncell)
{
    if (y >= 0 && y < ysize && x >= 0 && x < xsize) (*data)[x][y] = ncell;
}

void Lattice2D_no_inter::setF(int x, int y, int color, const array2D& nf)
{
    DistributionSetType2D f = (*data)[x][y].getF();
    f[color] = nf;
    (*data)[x][y].setF(f);
}

void Lattice2D_no_inter::setF(int x, int y, int color, int pos, double value)
{
    DistributionSetType2D f = (*data)[x][y].getF();
    f[color][pos] = value;
    (*data)[x][y].setF(f);
}

void Lattice2D_no_inter::setBoundaries(const Boundaries& newBound)
{
    bound = newBound;
    buildWalls();
}

void Lattice2D_no_inter::linearIndex(int index, int& x, int& y)const
{
    x = (index)%xsize;
    y = (index)/xsize;
}


//====================== BOUNDARY TREATMENT ======================

const bool Lattice2D_no_inter::isBoundary(int x, int y)const
{
    return  (x == 0 && bound.west.getType()!= periodic) || (x == xsize-1 && bound.east.getType()!= periodic) || (y == ysize-1 && bound.north.getType()!= periodic) || (y == 0 && bound.south.getType()!= periodic);
}

void Lattice2D_no_inter::buildWalls()
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

const Cell2D Lattice2D_no_inter::boundaryNorthPres(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::boundaryNorthVelo(const Cell2D& tmp, double uy)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(0,uy),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D_no_inter::boundarySouthPres(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::boundarySouthVelo(const Cell2D& tmp, double uy)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(0,uy),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D_no_inter::boundaryWestPres(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::boundaryWestVelo(const Cell2D& tmp, double ux)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(ux,0),Vector2D(0,0)}}, param.getPhi2D()));
}

const Cell2D Lattice2D_no_inter::boundaryEastPres(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::boundaryEastVelo(const Cell2D& tmp, double ux)const
{
    return Cell2D(eqDistro({{1.0,0}}, {{Vector2D(ux,0),Vector2D(0,0)}}, param.getPhi2D()));
}

// corners
const Cell2D Lattice2D_no_inter::cornerNorthWest(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::cornerNorthEast(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::cornerSouthWest(const Cell2D& tmp, ColSet rho)const
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

const Cell2D Lattice2D_no_inter::cornerSouthEast(const Cell2D& tmp, ColSet rho)const
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

Lattice2D_no_inter& Lattice2D_no_inter::operator=(const Lattice2D_no_inter& other)
{
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1]);
    this->setParams(other.getParams());
    this->setBoundaries(other.getBoundaries());
    this->setOffset(other.getOffset());

    return *this;
}

const bool Lattice2D_no_inter::operator==(const Lattice2D_no_inter& other)const
{
    bool exit = true;
    DimSet2D extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1])
    {
        ParamSet pOther = other.getParams();
        if (!(param == pOther)) exit = false;

        Boundaries bOther = other.getBoundaries();
        if (!(bound == bOther)) exit = false;

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

void Lattice2D_no_inter::streamAndBouncePull(Cell2D& tCell, const direction2D& dir)const
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

const array2D Lattice2D_no_inter::eqDistro_no_inter(double rho, const Vector2D& u)
{
    array2D feq;
    double usqr = u*u;
    for (int i=0; i<9; i++)
    {
        double scal = u*DIRECTION_2D[i];
        feq[i] = rho * WEIGHTS_2D[i] * ( 1 + 3 * scal + 4.5 * (scal*scal) - 1.5 * usqr);
    }
    return feq;
}
