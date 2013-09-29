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

void restart_file(const Lattice& l, const Preprocess& p, const Timetrack time, const string& filename){

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
    field data = l.getData();
    for (int x = 0; x<extent[0];x++){
        for (int y = 0; y<extent[1];y++){
            file.write(reinterpret_cast<char*> (&data[x][y]), sizeof(Cell));
        }
    }

    // write Timetrack
    int count = time.getCount();
    double factor = time.getFactor();
    double dtini = time.getDTini();
    vector<int> refinelist = time.getList();
    unsigned int vsize = refinelist.size();

    int maxCount = time.getMaxCount();
    double maxTime = time.getMaxTime();

    file.write(reinterpret_cast<char*> (&count), sizeof count);
    file.write(reinterpret_cast<char*> (&factor), sizeof factor);
    file.write(reinterpret_cast<char*> (&dtini), sizeof dtini);
    file.write(reinterpret_cast<char*> (&vsize), sizeof vsize);
    for(unsigned int i = 0; i< vsize; i++){
        int tmp = refinelist[i];
        file.write(reinterpret_cast<char*> (&tmp), sizeof(int));
    }
    file.write(reinterpret_cast<char*> (&maxCount), sizeof maxCount);
    file.write(reinterpret_cast<char*> (&maxTime), sizeof maxTime);

    double ReynoldsMax = p.getReynoldsMax();    
    double Morton = p.getMorton();
    double Eotvos = p.getEotvos();
    double resolution = p.getResolution();
    double rho_l = p.getRhoL();
    double gamma = p.getGamma();
    double diameter = p.getDiameter();
    double c_s = p.getSoundspeed();
    double sigma = p.getSigma();
    double g = p.getGPhys(); 

    file.write(reinterpret_cast<char*> (&ReynoldsMax), sizeof(double));
    file.write(reinterpret_cast<char*> (&Morton), sizeof(double));
    file.write(reinterpret_cast<char*> (&Eotvos), sizeof(double));
    file.write(reinterpret_cast<char*> (&resolution), sizeof(double));
    file.write(reinterpret_cast<char*> (&rho_l), sizeof(double));
    file.write(reinterpret_cast<char*> (&gamma), sizeof(double));
    file.write(reinterpret_cast<char*> (&diameter), sizeof(double));
    file.write(reinterpret_cast<char*> (&c_s), sizeof(double));

    file.write(reinterpret_cast<char*> (&sigma), sizeof(double));
    file.write(reinterpret_cast<char*> (&g), sizeof(double));

    file.close();
}

const bool restart_read(Lattice& outL, Preprocess& p, Timetrack& t, const string& filename)
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

        Cell tmpCell;
        field data(boost::extents[extent[0]][extent[1]]);
        for(int x = 0; x<extent[0];x++){
            for(int y=0;y<extent[1];y++){
                file.read((char*) &tmpCell, sizeof(Cell));
                data[x][y] = tmpCell;
            }
        }

        Timetrack time;

        int count;
        file.read((char*) &count, sizeof count);
        time.setCount(count);
        
        double factor;
        file.read((char*) &factor, sizeof factor);
        time.setFactor(factor);
        
        double dtini;
        file.read((char*) &dtini, sizeof dtini);
        time.setDTini(dtini);

        unsigned int vsize;
        file.read((char*) &vsize, sizeof vsize);
        
        vector<int> refinelist(vsize);
        for(unsigned int i = 0; i< vsize; i++){
            int tmp;
            file.read((char*) &tmp, sizeof(int));
            refinelist[i] = tmp; 
            }
        time.setVector(refinelist);

        int maxCount;
        file.read((char*) &maxCount, sizeof maxCount);
        time.setMaxCount(maxCount);
        double maxTime;
        file.read((char*) &maxTime, sizeof maxTime);
        time.setMaxTime(maxTime);


        double ReynoldsMax, Morton, Eotvos;
        double resolution, rho_l, gamma;
        double diameter, c_s, sigma,  g; 

        file.read((char*) &ReynoldsMax, sizeof(double));
        file.read((char*) &Morton, sizeof(double));
        file.read((char*) &Eotvos, sizeof(double));
        file.read((char*) &resolution, sizeof(double));
        file.read((char*) &rho_l, sizeof(double));
        file.read((char*) &gamma, sizeof(double));
        file.read((char*) &diameter, sizeof(double));
        file.read((char*) &c_s, sizeof(double));
        file.read((char*) &sigma, sizeof(double));
        file.read((char*) &g, sizeof(double));

        Preprocess prepro(ReynoldsMax, Morton, Eotvos, resolution, rho_l, gamma, diameter, c_s, sigma, g);
        
        file.close();
        outL.setParams(param);
        outL.setData(data, extent[0], extent[1]);
        outL.setTimetrack(time);
        t = time;
        p = prepro;
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
    paramLog << "speedlimit = "        << p.getSpeedlimit()<< " / m s^-1" << endl; 
    paramLog << "dt = "         << p.getDeltaT()    << " / m " << endl;
    paramLog << "gravity = "          << p.getG()         << " / -" << endl;

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
    // initialzing strings and fallback values
    map<string,double> mm;
        mm.insert(pair<string,double>("omega_red",1));
        mm.insert(pair<string,double>("omega_blue",1));
        mm.insert(pair<string,double>("rho_red",1));
        mm.insert(pair<string,double>("gamma",1000));
        mm.insert(pair<string,double>("alpha_blue",0.2));
        mm.insert(pair<string,double>("delta",0.1));
        mm.insert(pair<string,double>("beta",0.99));
        mm.insert(pair<string,double>("sigma",1e-4));
        mm.insert(pair<string,double>("speedlimit",1));
        mm.insert(pair<string,double>("timestep",0.001));
        mm.insert(pair<string,double>("g",9.81));


        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( inputQuery(filename,it->first,tmp) == true ) it->second = tmp;
        }

    ParamSet params(mm.at("omega_red"),mm.at("omega_blue"),mm.at("rho_red"),mm.at("gamma"),mm.at("sigma"),mm.at("g"), mm.at("speedlimit"),mm.at("timestep"),mm.at("alpha_blue"),mm.at("delta"),mm.at("beta"));
    return params;
}

const Preprocess getFilePreprocess(const string& filename){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
        mm.insert(pair<string,double>("Reynolds",50));
        mm.insert(pair<string,double>("Morton",1e-3));
        mm.insert(pair<string,double>("Eotvos",20));
        mm.insert(pair<string,double>("resolution",40));
        mm.insert(pair<string,double>("rho_l",1000));
        mm.insert(pair<string,double>("gamma",5));
        mm.insert(pair<string,double>("diameter",0.1));
        mm.insert(pair<string,double>("c_s",10));
        mm.insert(pair<string,double>("sigma",1e-4));
        mm.insert(pair<string,double>("g",10));


        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( inputQuery(filename,it->first,tmp) == true ) it->second = tmp;
        }

    Preprocess prepro(mm.at("Reynolds"),mm.at("Morton"),mm.at("Eotvos"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"),mm.at("diameter"),mm.at("c_s"),mm.at("sigma"), mm.at("g"));
    return prepro;
}

const Timetrack getFileTimetrack(const string& filename, const Preprocess& prepro){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("factor",1.1));
    mm.insert(pair<string,double>("max_steps",1e5));
    mm.insert(pair<string,double>("max_time",5));

    // cycling through the input file
    double tmp;
    for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
        if( inputQuery(filename,it->first,tmp) == true ) it->second = tmp;
    }

    Timetrack time(prepro.getTimestep(), mm.at("factor"), mm.at("max_steps"), mm.at("max_time"));
 
    return time;
}