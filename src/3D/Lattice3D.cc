#include"Lattice3D.h"

///////////////////////////// PUBLIC /////////////////////////////

//=========================== LIFECYCLE ===========================

Lattice3D::Lattice3D(int x_size, int y_size, int z_size, double fzero_dense, double fzero_dilute):
xsize(x_size)
,ysize(y_size)
,zsize(z_size)
,data(new field3D(boost::extents[xsize][ysize][zsize]))
,param()
,bound()
,bubblebox(0,0,0,xsize,ysize,zsize)
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
,bound(other.getBoundaries())
,bubblebox(other.getBubbleBox())
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
        // tmpDelta = GRAD_WEIGHTS_3D[q] * (*data)[ dir[q].x ][ dir[q].y ][ dir[q].z ].getDeltaRho();
        tmpDelta = GRAD_WEIGHTS_3D[q] * (*data)[ dir[q].x ][ dir[q].y ][ dir[q].z ].calcPsi();
        grad.x += DIRECTION_3D[q].x * tmpDelta;
        grad.y += DIRECTION_3D[q].y * tmpDelta;
        grad.z += DIRECTION_3D[q].z * tmpDelta;
    }
    return grad;
}

const Vector3D Lattice3D::getSurfaceNormal(int x, int y, int z) const
{
    const Vector3D grad = getGradient(x,y,z);
    const double abs = grad.Abs();
    Vector3D n(0,0,0);
    if (abs > 1e-3) n = grad * (1.0 / abs);
    return n;
}

void Lattice3D::streamAll(int threads)
{
    field3D *newData = new field3D(boost::extents[xsize][ysize][zsize]);

    omp_set_num_threads (threads);
    #pragma omp parallel
    {
        #pragma omp for collapse(3) schedule(dynamic, 100)       
        for (int z = 0; z < zsize; z++)
        {
            for (int y = 0; y < ysize; y++)
            {
                for (int x = 0; x < xsize; x++)
                {
                    const direction3D dir = directions(x,y,z);
                    Cell3D tmpCell = (*data)[x][y][z];
                   
                    if (tmpCell.getIsSolid() == false) streamAndBouncePull(tmpCell,dir);
                   
                    tmpCell.calcRho();
                    (*newData)[x][y][z] = tmpCell;
                }
            }
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
    const double dt = param.getDeltaT();

    double g(0);
    if(gravity == true) g = param.getG();

    const Matrix3D relaxation_matrix_outside(param.getRelaxation3D(1),false);
    const Matrix3D single_relax_outside = INV_TRAFO_MATRIX3D * relaxation_matrix_outside * TRAFO_MATRIX3D;

    #pragma omp parallel firstprivate(single_relax_outside)
    {
        #pragma omp for collapse(3) schedule(dynamic,100)
        for (int x=0; x<xsize; x++)
        {
            for (int y=0; y<ysize; y++)
            {
                for (int z=0; z<zsize; z++)
                {
                    Cell3D tmpCell = (*data)[x][y][z];
                    if(tmpCell.getIsSolid() == false && isBoundary(x,y,z) == false)
                    {

                        DistributionSetType3D  fTmp;
                        const DistributionSetType3D fCell = tmpCell.getF();

                        const ColSet rho_k = tmpCell.getRho();
                        const double rho = sum(rho_k);
                        VeloSet3D u = tmpCell.getU();
                        double fges;

                        if(bubblebox.isInBox(x,y,z) == false)
                        {
                            const DistributionSetType3D fEq = eqDistro(rho_k, u, phi);
                            const DistributionSetType3D diff = distro_diff_3D(fCell, fEq);
                            const DistributionSetType3D single_phase_col =  single_relax_outside * diff;
                            for (int q=0; q<19; q++)
                            {
                                for (int color=0;color<=1; color++)
                                {
                                    fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];
                                } // end for color
                                fges = fTmp[0][q]+fTmp[1][q];

                                if (rho > 0)
                                {
                                    fTmp[0][q] = rho_k[0]/rho  * fges;
                                    fTmp[1][q] = rho_k[1]/rho * fges;
                                } // end if rho
                            } // end for q
                        } // end if bubblebox
                        else 
                        {
                            const double psi = tmpCell.calcPsi();
                            const RelaxationPar3D relax = param.getRelaxation3D(psi);
                            const double omega = param.getOmega(psi);

                            const Matrix3D relaxation_matrix(relax,false);
                            const Matrix3D single_relax = INV_TRAFO_MATRIX3D * relaxation_matrix * TRAFO_MATRIX3D;
                            const Matrix3D forcing_factor(relax,true);    // (I - 0.5 S) -> ( 1 - 0.5 omega)
                            const Matrix3D forcing_relax = INV_TRAFO_MATRIX3D * forcing_factor * TRAFO_MATRIX3D;
                            const ColSet A_k = param.getAk(omega);

                            const Vector3D G(0 , 0, g*(rho - rho_k[0]));

                            if(gravity == true)
                            {
                                u[0] = u[0] + G *  (dt/(2* rho)) ;
                                u[1] = u[1] + G *  (dt/(2* rho)) ;
                            }

                            const DistributionSetType3D fEq = eqDistro(rho_k, u, phi);
                            const DistributionSetType3D diff = distro_diff_3D(fCell, fEq);
                            const DistributionSetType3D single_phase_col =  single_relax * diff;

                            const DistributionSetType3D second_forcing_term = forcing_relax * calculate_forcing_term(G,u); // M^{-1} F'
                            const Vector3D grad = getGradient(x,y,z);
                            const double av = grad.Abs();

                            double scal;
                            double recolor;

                            // q = 0
                            for (int color=0;color<=1; color++)
                            {
                                fTmp[color][0] =  fCell[color][0] - single_phase_col[color][0]; //single
                                if (gravity == true) fTmp[color][0] +=  dt * second_forcing_term[color][0]; // forcing term
                                fTmp[color][0] += A_k[color] * (- (av/2 * B_3D[0])); // perturbation
                            }
                            fges = fTmp[0][0]+fTmp[1][0];
                            if (rho > 0)
                            {
                                fTmp[0][0] = rho_k[0]/rho  * fges;
                                fTmp[1][0] = rho_k[1]/rho * fges;
                            }


                            for (int q=1; q<19; q++)
                            {
                                scal = grad*DIRECTION_3D[q];
                                double gradient_collision(0);
                                if (av > 0 ) gradient_collision = av/2 * (WEIGHTS_3D[q] * ( scal*scal )/(av*av) - B_3D[q]);

                                for (int color=0;color<=1; color++)
                                {
                                    fTmp[color][q] =  fCell[color][q] - single_phase_col[color][q];
                                    if (gravity == true) fTmp[color][q] +=  dt * second_forcing_term[color][q];
                                    fTmp[color][q] += A_k[color] * gradient_collision;
                                }

                                fges = fTmp[0][q]+fTmp[1][q];

                                // recoloring
                                if (rho > 0 && av>0) recolor = beta * (rho_k[0] * rho_k[1])/(rho*rho) * (scal /(av * DIRECTION_ABS_3D[q])) * (rho_k[0] * phi.at(0).at(q) + rho_k[1] * phi.at(1).at(q));
                                else recolor = 0;
                                if (rho > 0)
                                {
                                    fTmp[0][q] = rho_k[0]/rho  * fges + recolor;
                                    fTmp[1][q] = rho_k[1]/rho * fges - recolor;
                                }
                            } // end for
                        } // end else bubblebox

                        tmpCell.setF(fTmp);
                    } // end if solid
                    (*newData)[x][y][z] = tmpCell;                    
                } // end for z
            } // end for y
        } // end for x
    } // end pragma parallel

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

void Lattice3D::evaluateBoundaries(int threads)
{
    omp_set_num_threads (threads);

    #pragma omp parallel
    {
        // north boundary
        if(bound.north.getType() > 1)
        {
            int lowerX = 0;
            int upperX = xsize;

            int lowerY = 0;
            int upperY = ysize;

            if(bound.west.getType()  > 1) // north west corner
            {
                lowerY = 1;
            }
            if(bound.east.getType() > 1)
            {
                upperY = ysize - 1;
            }

            if(bound.back.getType()  > 1) // north west corner
            {
                lowerX = 1;
            }
            if(bound.front.getType() > 1)
            {
                upperX = xsize - 1;
            }

            if(bound.north.getType() == velocity)
            {
                double u_z = bound.north.getVelocity()[0].z;
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    for (int y = lowerY; y<upperY;y++)
                    {
                        (*data)[x][y][zsize-1] = boundaryVeloNS((*data)[x][y][zsize-1], u_z);
                    }
                }
            }
            if(bound.north.getType() == shear)
            {
                double u_0 = bound.west.getVelocity()[0].z;
                double u_1 = bound.east.getVelocity()[0].z;

                #pragma omp for schedule(static,10) nowait
                for (int y=lowerY; y<upperY; y++)
                {
                    double u_z = get_shearvelocity_x(y,u_0,u_1);
                    for (int x=lowerX; x<upperX; x++)
                    {
                        (*data)[x][y][zsize-1] = boundaryVeloNS((*data)[x][y][zsize-1], u_z);
                    }                    
                }
            }
        } // end north

        // south boundary
        if(bound.south.getType() > 1)
        {
            int lowerX = 0;
            int upperX = xsize;

            int lowerY = 0;
            int upperY = ysize;

    
            if(bound.west.getType()  > 1) // north west corner
            {
                lowerY = 1;
            }
            if(bound.east.getType() > 1)
            {
                upperY = ysize - 1;
            }

            if(bound.back.getType()  > 1) // north west corner
            {
                lowerX = 1;
            }
            if(bound.front.getType() > 1)
            {
                upperX = xsize - 1;
            }

            if(bound.north.getType() == velocity)
            {
                double u_z = bound.north.getVelocity()[0].z;
                #pragma omp for schedule(static,10) nowait
                for (int x=lowerX; x<upperX; x++)
                {
                    for (int y = lowerY; y<upperY;y++)
                    {
                        (*data)[x][y][0] = boundaryVeloNS((*data)[x][y][zsize-1], u_z);
                    }
                }
            }
            if(bound.north.getType() == shear)
            {
                double u_0 = bound.west.getVelocity()[0].z;
                double u_1 = bound.east.getVelocity()[0].z;

                #pragma omp for schedule(static,10) nowait
                for (int y=lowerY; y<upperY; y++)
                {
                    double u_z = get_shearvelocity_x(y,u_0,u_1);
                    for (int x=lowerX; x<upperX; x++)
                    {
                        (*data)[x][y][0] = boundaryVeloNS((*data)[x][y][zsize-1], u_z);
                    }                    
                }
            }
        } // end south 
    }
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

void Lattice3D::setShearProfile(double gradient, double offset)
{
    const double m = gradient;
    const double n = offset; // - (m * ( xsize *  param.getDeltaX() ) );
    for (int x = 0; x < xsize; x++)
    {
        VeloSet3D u;
        u[0] = Vector3D(0,m*x + n);
        u[1] = Vector3D(0,0);

        for (int y = 0; y < ysize; y++)
        {
            for (int z = 0; z < zsize; z++)
            {
                Cell3D tmpCell = (*data)[x][y][z];
                tmpCell.calcRho();
                const ColSet rho = tmpCell.getRho();
                const DistributionSetType3D phi = param.getPhi3D();

                tmpCell.setF(eqDistro(rho,u, phi));
                (*data)[x][y][z] = tmpCell;
            }
        }
    }
}

//=========================== ACCESSORS ===========================

const DimSet3D Lattice3D::getSize()const
{
    DimSet3D pony = {{xsize, ysize, zsize}}; 
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

void Lattice3D::setBoundaries(const Boundaries& newBound)
{
    bound = newBound;
    buildWalls();
}

void Lattice3D::linearIndex(int index, int& x, int& y, int& z)const
{
    x = (index)%xsize;
    y = ((index)/xsize)%ysize;
    z = (index)/(xsize * ysize);
}

//=========================== LATTICE CUTOUT ===========================

const std::vector<int> Lattice3D::findBubbleCells()const
{
    std::vector<int> indices;
    const int range = xsize * ysize * zsize;
    int x,y,z;

    for (int index = 0;  index < range; index++)
    {
        linearIndex(index,x,y,z);

        const Cell3D tmpCell = (*data)[x][y][z];

        if (tmpCell.calcPsi() < BUBBLE_CRITERION && tmpCell.getIsSolid() == false)
        {
            indices.push_back(index);
        }        
    }
    return indices;
}

void Lattice3D::copyCellsFromOther(const Lattice3D& other, const std::vector<int>& indices)
{
    int x,y,z;
    for(int index: indices)
    {
        linearIndex(index,x,y,z);
        (*data)[x][y][z] = other.getCell(x,y,z);
    }
}

const boost::array<Vector3D,2> Lattice3D::getBubbleData()const
{
    Vector3D momentum(0,0,0);
    Vector3D tmp_position(0,0,0);
    //Vector3D force(0,0,0);
    double rho_sum(0);


    for (int x = 0; x<xsize;x++)
    {
        for (int y = 0; y<ysize;y++)
        {
            for (int z = 0; z<zsize; z++)
            {
                const Cell3D tmp_cell = (*data)[x][y][z];
                const VeloSet3D tmp_velo = tmp_cell.getU();
                const ColSet tmp_rho = tmp_cell.getRho();
                if ( tmp_cell.calcPsi() < 0.99 && tmp_cell.getIsSolid() == false)
                {
                    momentum = momentum + (tmp_velo[1] * tmp_rho[1]);
                    tmp_position = tmp_position + (Vector3D(x,y,z) * tmp_rho[1]);
                    //force = force +  getSurfaceNormal(x, y, z) * getDivergence(x, y, z);
                    rho_sum += tmp_rho[1];
                }
            }
        }
    }

    const Vector3D velocity = momentum * (1.0/ rho_sum);
    const Vector3D position = tmp_position * (1.0/ rho_sum);
    
    // boost::array<Vector3D,3> result = {{position,velocity,force}};
    boost::array<Vector3D,2> result = {{position,velocity}};
    return result;

}

//====================== BOUNDARY TREATMENT ======================

const bool Lattice3D::isBoundary(int x, int y, int z)const
{
    return  (x == 0 && bound.back.getType()!= periodic) || (x == xsize-1 && bound.front.getType()!= periodic) || (y == 0 && bound.west.getType()!= periodic) || (y == ysize-1 && bound.east.getType()!= periodic) || (z == 0 && bound.south.getType()!= periodic) || (z == zsize-1 && bound.north.getType()!= periodic);
}

void Lattice3D::buildWalls()
{
    Cell3D wall(0,0,true);

    if(bound.south.getType() == bounceback)
    {
        Vector3D u3D = bound.south.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                (*data)[x][y][0] = wall;
            }
        }
    }

    if(bound.north.getType() == bounceback)
    {
        Vector3D u3D = bound.north.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            for (int y = 0; y < ysize; y++)
            {
                (*data)[x][y][zsize-1] = wall;
            }
        }
    }

    if(bound.west.getType() == bounceback)
    {
        Vector3D u3D = bound.west.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            for (int z = 0; z < zsize; z++)
            {
                (*data)[x][0][z] = wall;
            }
        }
    }

    if(bound.east.getType() == bounceback)
    {
        Vector3D u3D = bound.east.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int x=0; x<xsize; x++)
        {
            for (int z = 0; z < zsize; z++)
            {
                (*data)[x][ysize-1][z] = wall;
            }
        }
    }    

    if(bound.back.getType() == bounceback)
    {
        Vector3D u3D = bound.back.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            for (int z = 0; z < zsize; z++)
            {
                (*data)[0][y][z] = wall;
            }
        }
    }

    if(bound.front.getType() == bounceback)
    {
        Vector3D u3D = bound.front.getVelocity()[0];
        Vector3D u_wall(u3D.x,u3D.y,u3D.z);
        wall.setSolidVelocity(u_wall);

        for (int y=0; y<ysize; y++)
        {
            for (int z = 0; z < zsize; z++)
            {
                (*data)[xsize-1][y][z] = wall;
            }
        }
    }
}

const Cell3D Lattice3D::boundaryVeloNS(const Cell3D& tmp, double uz)const
{
    return Cell3D(eqDistro({{1.0,0}}, {{Vector3D(0,0,uz),Vector3D(0,0,0)}}, param.getPhi3D()));
}

const Cell3D Lattice3D::boundaryVeloEW(const Cell3D& tmp, double uy)const
{
    return Cell3D(eqDistro({{1.0,0}}, {{Vector3D(0,uy,0),Vector3D(0,0,0)}}, param.getPhi3D()));
}

const Cell3D Lattice3D::boundaryVeloFB(const Cell3D& tmp, double ux)const
{
    return Cell3D(eqDistro({{1.0,0}}, {{Vector3D(ux,0,0),Vector3D(0,0,0)}}, param.getPhi3D()));
}

//=========================== OPERATOR ===========================

Lattice3D& Lattice3D::operator=(const Lattice3D& other){
    this->setData(other.getData(), other.getSize()[0], other.getSize()[1], other.getSize()[2]);
    this->setParams(other.getParams());
    this->setBoundaries(other.getBoundaries());
    this->setBubbleBox(other.getBubbleBox());

    return *this;
}

const bool Lattice3D::operator==(const Lattice3D& other)const
{
    bool exit = true;
    DimSet3D extent = other.getSize();
    if (xsize == extent[0] && ysize == extent[1] && zsize == extent[2])
    {
        if (!(param == other.getParams())) exit = false;

        if (!(bound == other.getBoundaries())) exit = false;

        if (!(bubblebox == other.getBubbleBox())) exit = false;

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

void Lattice3D::streamAndBouncePull(Cell3D& tCell, const direction3D& dir)const
{
    const DistributionSetType3D f = tCell.getF();
    DistributionSetType3D ftmp;
    for (int color = 0; color<=1;color++)
    {
        ftmp[color][0] = (*data)[ dir[0].x ][ dir[0].y][ dir[0].z].getF()[color][0];

        for (int i=1;i<19;i++)
        {
            const Cell3D neighbor = (*data)[ dir[PULL_INDEX_3D[i]].x][ dir[PULL_INDEX_3D[i]].y][ dir[PULL_INDEX_3D[i]].z];
            if (neighbor.getIsSolid() == false)     // if(neighbor not solid?) -> stream
            {
                ftmp[color][i] = neighbor.getF()[color][i];
            }
            else // else -> bounce back
            {
                const ColSet rho = tCell.getRho();
                ftmp[color][i] = f[color][PULL_INDEX_3D[i]] - (2.0 * WEIGHTS_3D[i] * rho[color] * (DIRECTION_3D[i] * neighbor.getU()[0]) ) ; 
            } 
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