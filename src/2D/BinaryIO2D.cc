#include"BinaryIO2D.h"

//=========================== BINARY DUMP ===========================

void write_binary2D(const Lattice2D& l, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    ColSet extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    field2D data = l.getData();

    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell2D));
        }
    }
    file.close();
}

const bool read_binary2D(Lattice2D& outL, const string& filename){
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        ColSet extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Cell2D tmpCell;
        field2D data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell2D));
                data[x][y] = tmpCell;
            }
        }
        file.close();
        outL.setParams(param);
        outL.setData(data, extent[0], extent[1]);
    }
    else success = false;

    return success;
}

//=========================== RESTART FILES ===========================

void write_restart_file2D(const Lattice2D& l, const Preprocess& p, const Timetrack time, const string& filename){

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::out | ios::binary);
    file.seekp(0);

    // start to write
    ColSet extent = l.getSize();
    file.write(reinterpret_cast<char*> (&extent), sizeof extent);

    ParamSet param = l.getParams();
    file.write(reinterpret_cast<char*> (&param), sizeof param);

    // write the velocity distributions
    field2D data = l.getData();
    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell2D));
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
    double s_3 = p.getS_3();
    double s_5 = p.getS_5();
    double s_11 = p.getS_11();
    double s_17 = p.getS_17();
    int width = p.getWidth();
    int height = p.getHeight();

    file.write(reinterpret_cast<char*> (&ReynoldsMax), sizeof(double));
    file.write(reinterpret_cast<char*> (&Morton), sizeof(double));
    file.write(reinterpret_cast<char*> (&Eotvos), sizeof(double));
    file.write(reinterpret_cast<char*> (&resolution), sizeof(double));
    file.write(reinterpret_cast<char*> (&rho_l), sizeof(double));
    file.write(reinterpret_cast<char*> (&gamma), sizeof(double));
    file.write(reinterpret_cast<char*> (&mu_ratio), sizeof(double));    
    file.write(reinterpret_cast<char*> (&s_3), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_5), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_11), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_17), sizeof(double));
    file.write(reinterpret_cast<char*> (&width), sizeof(int));
    file.write(reinterpret_cast<char*> (&height), sizeof(int));

    file.close();
}

const bool read_restart_file2D(Lattice2D& outL, Preprocess& p, Timetrack& t, const string& filename)
{
    bool success;

    // setting up file
    fstream file(filename.c_str(),ios::in | ios::binary);
    if(file.is_open()){
        success = true;
        file.seekg(0);

        // start to read
        ColSet extent;
        file.read((char*) &extent, sizeof extent);

        ParamSet param;
        file.read((char*) &param, sizeof param);

        Cell2D tmpCell;
        field2D data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell2D));
                data[x][y] = tmpCell;
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
        double resolution, rho_l, gamma;
        double mu_ratio, s_3, s_5, s_11, s_17; 
        int width, height;

        file.read((char*) &ReynoldsMax, sizeof(double));
        file.read((char*) &Morton, sizeof(double));
        file.read((char*) &Eotvos, sizeof(double));
        file.read((char*) &resolution, sizeof(double));
        file.read((char*) &rho_l, sizeof(double));
        file.read((char*) &gamma, sizeof(double));
        file.read((char*) &mu_ratio, sizeof(double));
        file.read((char*) &s_3, sizeof(double));
        file.read((char*) &s_5, sizeof(double));
        file.read((char*) &s_11, sizeof(double));
        file.read((char*) &s_17, sizeof(double));

        file.read((char*) &width, sizeof(int));
        file.read((char*) &height, sizeof(int));

        Preprocess prepro(ReynoldsMax, Morton, Eotvos, resolution, rho_l, gamma, mu_ratio, s_3, s_5, s_11, s_17, width, height);
        
        file.close();
        outL.setParams(param);
        outL.setData(data, extent[0], extent[1]);
        t = time;
        p = prepro;
    }
    else success = false;

    return success;
}

//=========================== WRITE OUTPUT ===========================

void write_techplot_output2D(const Lattice2D& l, int iterNum)
{
    ofstream PsiFile;
    Cell2D tmp;
    ColSet extent = l.getSize();
    int xsize = static_cast<int> (extent[0]);
    int ysize = static_cast<int> (extent[1]);

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
    PsiFile << "\"ux2\""<<endl;
    PsiFile << "\"uy2\""<<endl;
    PsiFile << "\"u_abs\""<<endl;
    PsiFile << "\"grad_x\""<<endl;
    PsiFile << "\"grad_y\""<<endl;
    PsiFile << "ZONE T=\" "<< iterNum << "\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K=1,F=POINT"<< endl;

    PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
    
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            
            VeloSet2D u = tmp.getU();
            ColSet rho = tmp.getRho();
            Vector2D v;
            Vector2D gradient = l.getGradient(i, j);
            if(sum(rho) > 0) v = (u[0]*rho[0] + u[1]*rho[1]) * (1/sum(rho));

            PsiFile << i << "\t" << j << "\t" << "0 \t" << tmp.calcPsi() << "\t" << sum(rho) << "\t" << u[0].x << "\t" << u[0].y << "\t" << u[1].x << "\t" << u[1].y << "\t" << v.Abs() << "\t" << gradient.x << "\t" << gradient.y << endl;
        }
    }
    PsiFile.close();
}

void write_techplot_output_alternative2D(const Lattice2D& l, const string& filename)
{
    ofstream PsiFile;
    Cell2D tmp;
    ColSet extent = l.getSize();
    int xsize = static_cast<int> (extent[0]);
    int ysize = static_cast<int> (extent[1]);

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
    PsiFile << "\"ux2\""<<endl;
    PsiFile << "\"uy2\""<<endl;
    PsiFile << "\"u_abs\""<<endl;
    PsiFile << "\"grad_x\""<<endl;
    PsiFile << "\"grad_y\""<<endl;
    
    PsiFile << "ZONE T=\"0\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K=1,F=POINT"<< endl;
    PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
    
    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            double psi = tmp.calcPsi();
            PsiFile << i << "\t" << j << "\t" << "0 \t" << psi ;
            
            VeloSet2D u = tmp.getU();
            ColSet rho = tmp.getRho();
            Vector2D v;
            Vector2D gradient = l.getGradient(i, j);
            if(sum(rho) > 0) v = (u[0]*rho[0] + u[1]*rho[1]) * (1/sum(rho));
            if(psi < 1) {
                u[0].x = 0;
                u[0].y = 0;
            }
            if(psi > -0.99) {
                u[1].x = 0;
                u[1].y = 0;
            }
            PsiFile << "\t" << sum(rho) << "\t" << u[0].x*rho[0] << "\t" << u[0].y*rho[0] << "\t" << u[1].x*rho[1] << "\t" << u[1].y*rho[1] << "\t" << v.Abs() << "\t" << gradient.x << "\t" << gradient.y;
            PsiFile << endl;
        }
    }
    PsiFile.close();
}

void write_vtk_output2D(const Lattice2D& l, const string& filename)
{
    ofstream VTKFile;
    Cell2D tmp;
    int e;
    ColSet extent = l.getSize();
    int xsize = static_cast<int> (extent[0]);
    int ysize = static_cast<int> (extent[1]);

    // stringstream name;
    // name <<"test_"<< iterNum<<".vtk";

    // VTKFile.open( name.str().c_str());
    VTKFile.open(filename.c_str());

    VTKFile << "# vtk DataFile Version 3.1" << endl;
    VTKFile << "Lattice2D Boltzmann data" << endl;
    VTKFile << "ASCII" << endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << endl;

    VTKFile << "POINTS "<< (xsize+1) * (ysize+1)  <<" INT \n";

    for (int j = 0; j <= ysize; j++)
    {
        for (int i = 0; i <= xsize; i++)
        {
            VTKFile << i << " " << j << " 0 " ;
        }
        VTKFile<<endl;
    }

    VTKFile << "\nCELLS " << (xsize) * (ysize) << " " << (xsize) * (ysize) * 5 << "\n";
    
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            e = i+(xsize+1)*j;            
            VTKFile <<"4 "<< e << " " << e+1 << " "<< e + xsize +1 << " " << e + xsize + 2 << " ";
        }
        VTKFile<< endl;
    }
    VTKFile << "\nCELL_TYPES "<< (xsize) * (ysize) << "\n";
    for (int q = 0; q < (xsize * ysize); q++)
    {
        VTKFile <<"8 ";
    }

    VTKFile << "\nCELL_DATA "<< (xsize) * (ysize) << endl;

    VTKFile << "SCALARS Psi DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            VTKFile << tmp.calcPsi() << " ";
        }
    }

    VTKFile << "\nSCALARS Rho DOUBLE\nLOOKUP_TABLE default"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            VTKFile << sum(tmp.getRho()) << " ";
        }
    }

    VTKFile << "\nVECTORS j1 DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {        
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();

            VeloSet2D u = tmp.getU();
            ColSet rho = tmp.getRho();
            VTKFile << u[0].x * rho[0] << " " << u[0].y * rho[0] << " 0 ";
        }
    }

    VTKFile << "\nVECTORS j2 DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();

            VeloSet2D u = tmp.getU();
            ColSet rho = tmp.getRho();
            VTKFile << u[1].x * rho[1] << " " << u[1].y * rho[1] << " 0 ";
        }
    }

    VTKFile << "\nVECTORS gradient DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            Vector2D gradient = l.getGradient(i, j);    
            VTKFile << gradient.x << " " << gradient.y  << " 0 ";
        }
    }

    VTKFile.close();
}

