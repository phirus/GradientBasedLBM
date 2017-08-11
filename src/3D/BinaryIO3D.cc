#include"BinaryIO3D.h"

//=========================== BINARY DUMP ===========================

void write_binary3D(const Lattice3D& l, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    DimSet3D extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    Boundaries bound = l.getBoundaries();
    file.write(reinterpret_cast<char*> (&bound), sizeof bound);

    BubbleBox3D bubble = l.getBubbleBox();
    file.write(reinterpret_cast<char*> (&bubble), sizeof bubble);

    field3D data = l.getData();

    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            for (int z = 0; z<extent[2];z++)
            {
                file.write(reinterpret_cast<char*> (&data[x][y][z]), sizeof(Cell3D));
            }            
        }
    }
    file.close();
}

const bool read_binary3D(Lattice3D& outL, const string& filename){
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        DimSet3D extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Boundaries bound;
        file.read((char*) &bound, sizeof bound);

        BubbleBox3D bubble;
        file.read((char*) &bubble, sizeof bubble);

        Cell3D tmpCell;
        field3D data(boost::extents[extent[0]][extent[1]][extent[2]]);
        for(int x = 0; x<extent[0];x++)
        {
            for(int y=0;y<extent[1];y++)
            {
                for(int z=0;z<extent[2];z++)
                {
                    file.read((char*) &tmpCell, sizeof(Cell3D));
                    data[x][y][z] = tmpCell;
                }
            }
        }
        file.close();
        outL.setParams(param);
        outL.setBoundaries(bound);
        outL.setBubbleBox(bubble);
        outL.setData(data, extent[0], extent[1], extent[2]);
    }
    else success = false;

    return success;
}

//=========================== RESTART FILES ===========================

void write_restart_file3D(const Lattice3D& l, const Preprocess& p, const Timetrack time, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    DimSet3D extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    Boundaries bound = l.getBoundaries();
    file.write(reinterpret_cast<char*> (&bound), sizeof bound);

    BubbleBox3D bubble = l.getBubbleBox();
    file.write(reinterpret_cast<char*> (&bubble), sizeof bubble);

    // write the velocity distributions
    field3D data = l.getData();
    for (int x = 0; x<extent[0];x++)
    {
        for (int y = 0; y<extent[1];y++)
        {
            for (int z = 0; z<extent[2]; z++)
            {
                file.write(reinterpret_cast<char*> (&data[x][y][z]), sizeof(Cell3D));
            }
        }
    }

    // write Timetrack
    int count = time.getCount();
    int maxCount = time.getMaxCount();
    int output_interval = time.getOutputInt();
    int restart_interval = time.getRestartInt();

    file.write(reinterpret_cast<char*> (&count), sizeof count);
    file.write(reinterpret_cast<char*> (&maxCount), sizeof maxCount);
    file.write(reinterpret_cast<char*> (&output_interval), sizeof output_interval);
    file.write(reinterpret_cast<char*> (&restart_interval), sizeof restart_interval);

    // write Preprocess

    double ReynoldsMax = p.getReynoldsMax();    
    double Morton = p.getMorton();
    double Eotvos = p.getEotvos();
    double resolution = p.getResolution();
    double rho_l = p.getRhoL();
    double gamma = p.getGamma();
    double mu_ratio = p.getMuRatio();
    double bulk_visco = p.getBulkVisco();
    double s_3 = p.getS_3();
    double s_5 = p.getS_5();
    double s_11 = p.getS_11();
    double s_17 = p.getS_17();
    bool isShearFlow = p.getIsShearFlow();
    double shearRate = p.getShearRate();
    int xCells = p.getXCells();
    int yCells = p.getYCells();
    int zCells = p.getZCells();

    file.write(reinterpret_cast<char*> (&ReynoldsMax), sizeof(double));
    file.write(reinterpret_cast<char*> (&Morton), sizeof(double));
    file.write(reinterpret_cast<char*> (&Eotvos), sizeof(double));
    file.write(reinterpret_cast<char*> (&resolution), sizeof(double));
    file.write(reinterpret_cast<char*> (&rho_l), sizeof(double));
    file.write(reinterpret_cast<char*> (&gamma), sizeof(double));
    file.write(reinterpret_cast<char*> (&mu_ratio), sizeof(double));
    file.write(reinterpret_cast<char*> (&bulk_visco), sizeof(double));    
    file.write(reinterpret_cast<char*> (&s_3), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_5), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_11), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_17), sizeof(double));
    file.write(reinterpret_cast<char*> (&isShearFlow), sizeof(bool));
    file.write(reinterpret_cast<char*> (&shearRate), sizeof(double));
    file.write(reinterpret_cast<char*> (&xCells), sizeof(int));
    file.write(reinterpret_cast<char*> (&yCells), sizeof(int));
    file.write(reinterpret_cast<char*> (&zCells), sizeof(int));

    file.close();
}

const bool read_restart_file3D(Lattice3D& outL, Preprocess& p, Timetrack& t, const string& filename)
{
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        DimSet3D extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Boundaries bound;
        file.read((char*) &bound, sizeof bound);

        BubbleBox3D bubble;
        file.read((char*) &bubble, sizeof bubble);

        Cell3D tmpCell;
        field3D data(boost::extents[extent[0]][extent[1]][extent[2]]);
        for(int x = 0; x<extent[0];x++)
        {
            for(int y=0;y<extent[1];y++)
            {
                for(int z=0;z<extent[2];z++)
                {
                    file.read((char*) &tmpCell, sizeof(Cell3D));
                    data[x][y][z] = tmpCell;
                }
            }
        }

        int count;
        int maxCount;
        int output_interval;
        int restart_interval;

        file.read((char*) &count, sizeof count);
        file.read((char*) &maxCount, sizeof maxCount);
        file.read((char*) &output_interval, sizeof output_interval);
        file.read((char*) &restart_interval, sizeof restart_interval);

        Timetrack time(maxCount, output_interval, restart_interval);
        time.setCount(count);

        double ReynoldsMax, Morton, Eotvos;
        double resolution, rho_l, gamma, mu_ratio;
        double bulk_visco, s_3, s_5, s_11, s_17;
        bool isShearFlow;
        double shearRate; 
        int xCells, yCells, zCells;

        file.read((char*) &ReynoldsMax, sizeof(double));
        file.read((char*) &Morton, sizeof(double));
        file.read((char*) &Eotvos, sizeof(double));
        file.read((char*) &resolution, sizeof(double));
        file.read((char*) &rho_l, sizeof(double));
        file.read((char*) &gamma, sizeof(double));
        file.read((char*) &mu_ratio, sizeof(double));
        file.read((char*) &bulk_visco, sizeof(double));
        file.read((char*) &s_3, sizeof(double));
        file.read((char*) &s_5, sizeof(double));
        file.read((char*) &s_11, sizeof(double));
        file.read((char*) &s_17, sizeof(double));
        file.read((char*) &isShearFlow, sizeof(bool));
        file.read((char*) &shearRate, sizeof(double));
        file.read((char*) &xCells, sizeof(int));
        file.read((char*) &yCells, sizeof(int));
        file.read((char*) &zCells, sizeof(int));

        Preprocess prepro(ReynoldsMax, Morton, Eotvos, resolution, rho_l, gamma, mu_ratio, bulk_visco, s_3, s_5, s_11, s_17, isShearFlow, shearRate, xCells, yCells, zCells);
        
        file.close();
        outL.setParams(param);
        outL.setBoundaries(bound);
        outL.setBubbleBox(bubble);
        outL.setData(data, extent[0], extent[1], extent[2]);
        t = time;
        p = prepro;
    }
    else success = false;

    return success;
}

//=========================== WRITE OUTPUT ===========================

void write_techplot_output3D(const Lattice3D& l, int iterNum)
{
    ofstream PsiFile;
    Cell3D tmp;
    DimSet3D extent = l.getSize();
    int xsize = extent[0];
    int ysize = extent[1];
    int zsize = extent[2];

    stringstream name;
    name <<"psi_"<< iterNum<<".dat";

    PsiFile.open( name.str().c_str() );
    PsiFile << "TITLE = \" NewDataSet \" "<<endl;
    PsiFile << "Variables = \" x \" "<<endl;
    PsiFile << "\"y\""<<endl;
    PsiFile << "\"z\""<<endl;
    PsiFile << "\"psi\""<<endl;
    PsiFile << "\"rho\""<<endl;
    PsiFile << "\"ux1\""<<endl;
    PsiFile << "\"uy1\""<<endl;
    PsiFile << "\"uz1\""<<endl;
    PsiFile << "\"ux2\""<<endl;
    PsiFile << "\"uy2\""<<endl;
    PsiFile << "\"uz2\""<<endl;
    PsiFile << "\"u_abs\""<<endl;
    PsiFile << "\"grad_x\""<<endl;
    PsiFile << "\"grad_y\""<<endl;
    PsiFile << "\"grad_z\""<<endl;
    PsiFile << "ZONE T=\" "<< iterNum << "\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K="<<zsize<<", F=POINT"<< endl;

    PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
    
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            for (int k=0; k<zsize; k++)
            {
                tmp = l.getCell(i,j,k);
                tmp.calcRho();

                VeloSet3D u = tmp.getU();
                ColSet rho = tmp.getRho();
                Vector3D v;
                Vector3D gradient = l.getGradient(i, j, k);
                if(sum(rho) > 0) v = (u[0]*rho[0] + u[1]*rho[1]) * (1/sum(rho));

                PsiFile << i << "\t" << j << "\t" << k <<"\t" << tmp.calcPsi() << "\t" << sum(rho) << "\t" << u[0].x << "\t" << u[0].y << "\t"<< u[0].z << "\t" << u[1].x << "\t" << u[1].y << "\t"<< u[1].z << "\t" << v.Abs() << "\t" << gradient.x << "\t" << gradient.y << "\t" << gradient.z << endl;
            }
         }
    }
    PsiFile.close();
}

void write_techplot_output_alternative3D(const Lattice3D& l, const string& filename)
{
    ofstream PsiFile;
    Cell3D tmp;
    DimSet3D extent = l.getSize();
    int xsize = extent[0];
    int ysize = extent[1];
    int zsize = extent[2];

    stringstream name;
    name <<filename;

    PsiFile.open( name.str().c_str() );
    PsiFile << "TITLE = \" NewDataSet \" "<<endl;
    PsiFile << "Variables = \" x \" "<<endl;
    PsiFile << "\"y\""<<endl;
    PsiFile << "\"z\""<<endl;
    PsiFile << "\"psi\""<<endl;
    PsiFile << "\"rho\""<<endl;
    PsiFile << "\"ux1\""<<endl;
    PsiFile << "\"uy1\""<<endl;
    PsiFile << "\"uz1\""<<endl;
    PsiFile << "\"ux2\""<<endl;
    PsiFile << "\"uy2\""<<endl;
    PsiFile << "\"uz2\""<<endl;
    PsiFile << "\"u_abs\""<<endl;
    PsiFile << "\"grad_x\""<<endl;
    PsiFile << "\"grad_y\""<<endl;
    PsiFile << "\"grad_z\""<<endl;
    
    PsiFile << "ZONE T=\"0\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K=" <<zsize<<",F=POINT"<< endl;
    PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
    
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            for (int k=0; k<zsize; k++)
            {
                tmp = l.getCell(i,j,k);
                tmp.calcRho();
                double psi = tmp.calcPsi();
                PsiFile << i << "\t" << j << "\t" << k << "\t" << psi ;

                VeloSet3D u = tmp.getU();
                ColSet rho = tmp.getRho();
                Vector3D v;
                Vector3D gradient = l.getGradient(i, j, k);
                if(sum(rho) > 0) v = (u[0]*rho[0] + u[1]*rho[1]) * (1/sum(rho));

                if(psi < 1) 
                {
                    u[0].x = 0;
                    u[0].y = 0;
                    u[0].z = 0;
                }

                if(psi > -0.99) 
                {
                    u[1].x = 0;
                    u[1].y = 0;
                    u[1].z = 0;
                }

                PsiFile << "\t" << sum(rho) << "\t" << u[0].x*rho[0] << "\t" << u[0].y*rho[0] << "\t" << u[0].z*rho[0] << "\t" << u[1].x*rho[1] << "\t" << u[1].y*rho[1]  << "\t" << u[1].z*rho[1] << "\t" << v.Abs() << "\t" << gradient.x << "\t" << gradient.y << "\t" << gradient.z;
                PsiFile << endl;
            }
        }
    }
    PsiFile.close();
}

void write_vtk_output3D(const Lattice3D& l, const string& filename)
{
    ofstream VTKFile;
    Cell3D tmp;
    DimSet3D extent = l.getSize();
    int xsize = extent[0];
    int ysize = extent[1];
    int zsize = extent[2];

    VTKFile.open(filename.c_str());

    VTKFile << "# vtk DataFile Version 3.1" << endl;
    VTKFile << "Lattice Boltzmann data" << endl;
    VTKFile << "ASCII" << endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << endl;

    VTKFile << "POINTS "<< (xsize+1) * (ysize+1) * (zsize+1)  <<" INT \n";

    for (int k = 0; k<=zsize; k++)
    {
        for (int j = 0; j <= ysize; j++)
        {
            for (int i = 0; i <= xsize; i++)
            {
                VTKFile << i << " " << j << " " << k << "  " ;
            }
        }
        VTKFile<<endl;
    }

    VTKFile << "\nCELLS " << (xsize) * (ysize) * (zsize) << " " << (xsize) * (ysize) * (zsize) * 9 << "\n";
    
    int x,y,z;
    int e,f;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        e = x+ (xsize+1)*y + (xsize+1)*(ysize+1)*z;
        f = (xsize+1)*(ysize+1);
        VTKFile <<"8 "<< e << " " << e+1 << " "<< e + xsize +1 << " " << e + xsize + 2 << " " << e + f << " " << e + f + 1 << " "<< e + f + xsize +1 << " " << e + f + xsize + 2 << " ";
    }

    VTKFile<< endl;

    VTKFile << "\nCELL_TYPES "<< (xsize * ysize * zsize) << "\n";
    for (int q = 0; q < (xsize * ysize* zsize); q++)
    {
        VTKFile <<"11 ";
    }

    VTKFile << "\nCELL_DATA "<< (xsize * ysize * zsize) << endl;

    VTKFile << "SCALARS Psi DOUBLE\nLOOKUP_TABLE default"<<endl;

    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
        VTKFile << tmp.calcPsi() << " ";
    }

    VTKFile << "\nSCALARS Rho DOUBLE\nLOOKUP_TABLE default"<<endl;

    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
           VTKFile << sum(tmp.getRho()) << " ";
    }

    VTKFile << "\nVECTORS j1 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        ColSet rho = tmp.getRho();
        VTKFile << u[0].x * rho[0] << " " << u[0].y * rho[0] << " " << u[0].z * rho[0] << " ";
    }

    VTKFile << "\nVECTORS j2 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        ColSet rho = tmp.getRho();
        VTKFile << u[1].x * rho[1] << " " << u[1].y * rho[1] << " " << u[1].z * rho[1] << " ";
    }

    VTKFile << "\nVECTORS gradient DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        Vector3D gradient = l.getGradient(x, y, z);
        VTKFile << gradient.x << " " << gradient.y  << " " << gradient.z << " ";
    }

    VTKFile << "\nVECTORS u1 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        VTKFile << u[0].x << " " << u[0].y << " " << u[0].z << " ";
    }

    VTKFile << "\nVECTORS u2 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        VTKFile << u[1].x << " " << u[1].y << " " << u[1].z << " ";
    }

    VTKFile.close();
}

void write_vtk_output3D_verboose(const Lattice3D& l, const string& filename)
{
    ofstream VTKFile;
    Cell3D tmp;
    DimSet3D extent = l.getSize();
    int xsize = extent[0];
    int ysize = extent[1];
    int zsize = extent[2];

    VTKFile.open(filename.c_str());

    VTKFile << "# vtk DataFile Version 3.1" << endl;
    VTKFile << "Lattice Boltzmann data" << endl;
    VTKFile << "ASCII" << endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << endl;

    VTKFile << "POINTS "<< (xsize+1) * (ysize+1) * (zsize+1)  <<" INT \n";

    for (int k = 0; k<=zsize; k++)
    {
        for (int j = 0; j <= ysize; j++)
        {
            for (int i = 0; i <= xsize; i++)
            {
                VTKFile << i << " " << j << " " << k << "  " ;
            }
        }
        VTKFile<<endl;
    }

    VTKFile << "\nCELLS " << (xsize) * (ysize) * (zsize) << " " << (xsize) * (ysize) * (zsize) * 9 << "\n";
    
    int x,y,z;
    int e,f;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        e = x+ (xsize+1)*y + (xsize+1)*(ysize+1)*z;
        f = (xsize+1)*(ysize+1);
        VTKFile <<"8 "<< e << " " << e+1 << " "<< e + xsize +1 << " " << e + xsize + 2 << " " << e + f << " " << e + f + 1 << " "<< e + f + xsize +1 << " " << e + f + xsize + 2 << " ";
    }

    VTKFile<< endl;

    VTKFile << "\nCELL_TYPES "<< (xsize * ysize * zsize) << "\n";
    for (int q = 0; q < (xsize * ysize* zsize); q++)
    {
        VTKFile <<"11 ";
    }

    VTKFile << "\nCELL_DATA "<< (xsize * ysize * zsize) << endl;

    VTKFile << "SCALARS Psi DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
        VTKFile << tmp.calcPsi() << " ";
    }

    VTKFile << "\nSCALARS Rho DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
        VTKFile << sum(tmp.getRho()) << " ";
    }

    VTKFile << "\nSCALARS Rho1 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
        VTKFile << tmp.getRho()[0] << " ";
    }

    VTKFile << "\nSCALARS Rho2 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();
        VTKFile << tmp.getRho()[1] << " ";
    }

    // VTKFile << "\nTENSORS p1 DOUBLE"<<endl;
    // for (int index = 0; index < (xsize*ysize*zsize); index++)
    // {
    //     l.linearIndex(index, x, y, z);
    //     tmp = l.getCell(x,y,z);
    //     const pressureTensor3D p1 = tmp.getPressureTensor(0);
    //     VTKFile << p1.xx << " " << p1.xy << " " << p1.xz << endl;
    //     VTKFile << p1.yx << " " << p1.yy << " " << p1.yz << endl;
    //     VTKFile << p1.zx << " " << p1.zy << " " << p1.zz << endl;
    //     VTKFile << endl;
    // }

    // VTKFile << "\nTENSORS p2 DOUBLE"<<endl;
    // for (int index = 0; index < (xsize*ysize*zsize); index++)
    // {
    //     l.linearIndex(index, x, y, z);
    //     tmp = l.getCell(x,y,z);
    //     const pressureTensor3D p2 = tmp.getPressureTensor(1);
    //     VTKFile << p2.xx << " " << p2.xy << " " << p2.xz << endl;
    //     VTKFile << p2.yx << " " << p2.yy << " " << p2.yz << endl;
    //     VTKFile << p2.zx << " " << p2.zy << " " << p2.zz << endl;
    //     VTKFile << endl;
    // }

    VTKFile << "\nSCALARS p_trace1 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        const pressureTensor3D p = tmp.getPressureTensor(0);

        VTKFile << p.getTrace()/3.0 << " ";
    }

    VTKFile << "\nSCALARS p_trace2 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        const pressureTensor3D p = tmp.getPressureTensor(1);

        VTKFile << p.getTrace()/3.0 << " ";
    }

    VTKFile << "\nSCALARS p_1_correction DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        double rho = tmp.getRho()[0];
        Vector3D u = tmp.getU()[0];

        VTKFile << rho*(u*u) << " ";
    }

    VTKFile << "\nSCALARS p_2_correction DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        double rho = tmp.getRho()[1];
        Vector3D u = tmp.getU()[1];

        VTKFile << rho*(u*u) << " ";
    }

    VTKFile << "\nSCALARS p_1 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        double rho = tmp.getRho()[0];
        Vector3D u = tmp.getU()[0];
        const pressureTensor3D p = tmp.getPressureTensor(0);

        VTKFile << (p.getTrace()/3.0) - rho*(u*u) << " ";
    }

    VTKFile << "\nSCALARS p_2 DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        double rho = tmp.getRho()[1];
        Vector3D u = tmp.getU()[1];
        const pressureTensor3D p = tmp.getPressureTensor(1);

        VTKFile << (p.getTrace()/3.0) - rho*(u*u) << " ";
    }

    VTKFile << "\nSCALARS p_total DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        double rho1 = tmp.getRho()[0];
        double rho2 = tmp.getRho()[1];
        Vector3D u1 = tmp.getU()[0];
        Vector3D u2 = tmp.getU()[1];
        const pressureTensor3D p1 = tmp.getPressureTensor(0);
        const pressureTensor3D p2 = tmp.getPressureTensor(1);

        VTKFile << (p1.getTrace()/3.0) - rho1*(u1*u1) +  (p2.getTrace()/3.0) - rho2*(u2*u2) << " ";
    }

    VTKFile << "\nVECTORS j1 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        ColSet rho = tmp.getRho();
        VTKFile << u[0].x * rho[0] << " " << u[0].y * rho[0] << " " << u[0].z * rho[0] << " ";
    }

    VTKFile << "\nVECTORS j2 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        ColSet rho = tmp.getRho();
        VTKFile << u[1].x * rho[1] << " " << u[1].y * rho[1] << " " << u[1].z * rho[1] << " ";
    }

    VTKFile << "\nVECTORS gradient DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        Vector3D gradient = l.getGradient(x, y, z);
        VTKFile << gradient.x << " " << gradient.y  << " " << gradient.z << " ";
    }

    VTKFile << "\nVECTORS u1 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        VTKFile << u[0].x << " " << u[0].y << " " << u[0].z << " ";
    }

    VTKFile << "\nVECTORS u2 DOUBLE"<<endl;
    for (int index = 0; index < (xsize*ysize*zsize); index++)
    {
        l.linearIndex(index, x, y, z);
        tmp = l.getCell(x,y,z);
        tmp.calcRho();

        VeloSet3D u = tmp.getU();
        VTKFile << u[1].x << " " << u[1].y << " " << u[1].z << " ";
    }

    VTKFile.close();
}

void writeBubbleFitData(const Lattice3D& l, const string& filename)
{
    ofstream BubbleFit;
    BubbleFit.open( filename.c_str() );
    BubbleFit << "x;y;z;rho_b\n";

    std::vector<int> indices = l.findBubbleCells();

    for(int index : indices)
    {
        int x,y,z;
        l.linearIndex(index, x, y, z);
        Cell3D tmp = l.getCell(x, y, z);

        BubbleFit << x << ";" << y << ";" << z << ";" << tmp.getRho()[1] << "\n";
    }
    BubbleFit.close();
}
