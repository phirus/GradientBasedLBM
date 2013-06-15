#include"binaryIO.h"

void binary_output(const Lattice& l, const string& filename){

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

    field data = l.getData();

    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell));
        }
    }

    file.close();
}

const bool binary_input(Lattice& outL, const string& filename){
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

        Cell tmpCell;
        field data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell));
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

void techplotOutput(const Lattice& l, int iterNum, bool verbose)
{
    ofstream PsiFile;
    Cell tmp;
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
    if (verbose == true)
    {
        PsiFile << "\"rho\""<<endl;
        PsiFile << "\"ux\""<<endl;
        PsiFile << "\"uy\""<<endl;
        PsiFile << "\"u_abs\""<<endl;
    }
    PsiFile << "ZONE T=\" "<< iterNum << "\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K=1,F=POINT"<< endl;
    if (verbose == true) PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
    else PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;

    for (int j=0; j<ysize; j++)
    {
        for (int i=0; i<xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            PsiFile << i << "\t" << j << "\t" << "0 \t" << tmp.calcPsi() ;
            if (verbose == true)
            {
                Vector u = tmp.getU();
                PsiFile << "\t" << sum(tmp.getRho()) << "\t" << u.x << "\t" << u.y << "\t" << u.abs() ;
            }
            PsiFile << endl;
        }
    }
    PsiFile.close();
}

void vtkOutput(const Lattice& l, int iterNum)
{
    ofstream VTKFile;
    Cell tmp;
    int e;
    ColSet extent = l.getSize();
    int xsize = static_cast<int> (extent[0]);
    int ysize = static_cast<int> (extent[1]);

    stringstream name;
    name <<"test_"<< iterNum<<".vtk";

    VTKFile.open( name.str().c_str() );
    VTKFile << "# vtk DataFile Version 3.1" << endl;
    VTKFile << "Lattice Boltzmann data" << endl;
    VTKFile << "ASCII" << endl;
    VTKFile << "DATASET UNSTRUCTURED_GRID" << endl;

    VTKFile << "POINTS "<< (xsize+1) * (ysize+1)  <<" INT "<< endl;

    for (int j = 0; j <= ysize; j++)
    {
        for (int i = 0; i <= xsize; i++)
        {
            VTKFile << i << " " << j << " 0 " ;
        }
        VTKFile<<endl;
    }

    VTKFile << "\nCELLS " << (xsize) * (ysize) << " " << (xsize) * (ysize) * 5 << endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            e = i+(xsize+1)*j;
            VTKFile <<"4 "<< e << " " << e+1 << " "<< e + xsize + 2 << " "<< e + xsize +1 << " ";
        }
        VTKFile<< endl;
    }
    VTKFile << "\nCELL_TYPES "<< (xsize) * (ysize) << endl;
    for (int q = 0; q < (xsize * ysize); q++)
    {
        VTKFile <<"9 ";
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

    VTKFile << "\nVECTORS u DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();
            Vector u = tmp.getU();
            VTKFile << u.x << " " << u.y << " 0 ";
        }
    }
    VTKFile.close();
}

void paramLogOut(const Lattice& l){
    ofstream paramLog;
    ColSet extent = l.getSize();
    ParamSet p = l.getParams();
    int xsize = static_cast<int> (extent[0]);
    int ysize = static_cast<int> (extent[1]);

    stringstream name;
    name <<"paramLog";

    paramLog.open( name.str().c_str() );
    paramLog << "# used setup parameters\n" << endl;
    paramLog << "xsize = "      << xsize            << " Cells " << endl;
    paramLog << "ysize = "      << ysize            << " Cells " << endl;
    paramLog << "omega_red = "  << p.getOmegaRed()  << " /-" << endl;
    paramLog << "omega_blue = " << p.getOmegaBlue() << " /-" << endl;
    paramLog << "rho_red = "    << p.getRhoR()      << " / kg m^-3" << endl;
    paramLog << "gamma = "      << p.getGamma()     << " /-" << endl;
    paramLog << "alpha_blue = " << p.getAlpha()     << " /-" << endl;
    paramLog << "delta = "      << p.getInterfaceThickness() << " /-" << endl;
    paramLog << "beta = "       << p.getBeta()      << " /-" << endl;
    paramLog << "sigma = "      << p.getSigma()     << " /?" << endl; //TODO get SI units of sigma
    paramLog << "c_s = "        << p.getSoundspeed()<< " / m s^-1" << endl; 
    paramLog << "dx = "         << p.getDeltaX()    << " / m " << endl;

    paramLog << "\n# successive parameters" << endl;
    paramLog << "dt = "         << p.getDeltaT()    << " / s" << endl;
    paramLog << "g = "          << p.getG()         << " /?" << endl;



    paramLog.close();
}

const bool inputQuery(const string& filename, const string& query, double& value){
    bool success = false;
    value = 0;
    string lineString;
    string word;
    fstream file(filename.c_str(),ios::in);

    while(file.good()){
    getline(file,lineString);
    if ( lineString.substr(0,1) != "#" ){
        istringstream lineStream(lineString);
        lineStream >> word;
        if (word == query) {
            lineStream >> word; // gets the "="
            lineStream >> value;
            success = true;
            break;
            }         
        }
    }
    file.close();
    return success;
}

const ParamSet getFileParams(const string& filename){
    vector<string> tags;
    vector<double> val;
    
    tags.push_back("omega_red");
    val.push_back(1);

    tags.push_back("omega_blue");
    val.push_back(1);
    
    tags.push_back("rho_red"); 
    val.push_back(1);
    
    tags.push_back("gamma");
    val.push_back(1000);
    
    tags.push_back("alpha_blue");
    val.push_back(0.2);
    
    tags.push_back("delta");
    val.push_back(0.1);
    
    tags.push_back("beta");
    val.push_back(0.99);
    
    tags.push_back("sigma");
    val.push_back(1e-4);
    
    tags.push_back("c_s");
    val.push_back(1484);
    
    tags.push_back("dx"); 
    val.push_back(0.001);

    double tmp;
    int imax = tags.size();
    for(int i = 0; i<imax; i++){
        if( inputQuery(filename,tags[i],tmp) == true ) val[i] = tmp;
    }

    ParamSet params(val[0],val[1],val[2],val[3],val[4],val[5],val[6],val[7],val[8],val[9]);
    return params;
}