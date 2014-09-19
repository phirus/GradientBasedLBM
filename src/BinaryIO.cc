#include"BinaryIO.h"

//=========================== BINARY DUMP ===========================

void write_binary(const Lattice& l, const string& filename){

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

const bool read_binary(Lattice& outL, const string& filename){
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

//=========================== RESTART FILES ===========================

void write_restart_file(const Lattice& l, const Preprocess& p, const Timetrack time, const string& filename){

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
    // double maxTime = time.getMaxTime();
    int techplot_interval = time.getTechPlotInt();
    int restart_interval = time.getRestartInt();

    file.write(reinterpret_cast<char*> (&count), sizeof count);
    file.write(reinterpret_cast<char*> (&factor), sizeof factor);
    file.write(reinterpret_cast<char*> (&dtini), sizeof dtini);
    file.write(reinterpret_cast<char*> (&vsize), sizeof vsize);
    for(unsigned int i = 0; i< vsize; i++){
        int tmp = refinelist[i];
        file.write(reinterpret_cast<char*> (&tmp), sizeof(int));
    }
    file.write(reinterpret_cast<char*> (&maxCount), sizeof maxCount);
    // file.write(reinterpret_cast<char*> (&maxTime), sizeof maxTime);
    file.write(reinterpret_cast<char*> (&techplot_interval), sizeof techplot_interval);
    file.write(reinterpret_cast<char*> (&restart_interval), sizeof restart_interval);

    double ReynoldsMax = p.getReynoldsMax();    
    double Morton = p.getMorton();
    double Eotvos = p.getEotvos();
    double resolution = p.getResolution();
    double rho_l = p.getRhoL();
    double gamma = p.getGamma();
    double mu_ratio = p.getMuRatio();
    double s_3 = p.getS_3();
    double s_5 = p.getS_5();
    int width = p.getWidth();
    int height = p.getHeight();

    file.write(reinterpret_cast<char*> (&ReynoldsMax), sizeof(double));
    file.write(reinterpret_cast<char*> (&Morton), sizeof(double));
    file.write(reinterpret_cast<char*> (&Eotvos), sizeof(double));
    file.write(reinterpret_cast<char*> (&resolution), sizeof(double));
    file.write(reinterpret_cast<char*> (&rho_l), sizeof(double));
    file.write(reinterpret_cast<char*> (&gamma), sizeof(double));
    // file.write(reinterpret_cast<char*> (&diameter), sizeof(double));
    file.write(reinterpret_cast<char*> (&mu_ratio), sizeof(double));    
    // file.write(reinterpret_cast<char*> (&c_s), sizeof(double));
    // file.write(reinterpret_cast<char*> (&sigma), sizeof(double));
    // file.write(reinterpret_cast<char*> (&g), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_3), sizeof(double));
    file.write(reinterpret_cast<char*> (&s_5), sizeof(double));
    file.write(reinterpret_cast<char*> (&width), sizeof(int));
    file.write(reinterpret_cast<char*> (&height), sizeof(int));

    file.close();
}

const bool read_restart_file(Lattice& outL, Preprocess& p, Timetrack& t, const string& filename)
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

        // double maxTime;
        // file.read((char*) &maxTime, sizeof maxTime);
        // time.setMaxTime(maxTime);

        int techplot_interval;
        file.read((char*) &techplot_interval, sizeof techplot_interval);
        time.setTechPlotInt(techplot_interval);

        int restart_interval;
        file.read((char*) &restart_interval, sizeof restart_interval);
        time.setRestartInt(restart_interval);

        double ReynoldsMax, Morton, Eotvos;
        double resolution, rho_l, gamma;
        double mu_ratio, s_3, s_5; 
        int width, height;

        file.read((char*) &ReynoldsMax, sizeof(double));
        file.read((char*) &Morton, sizeof(double));
        file.read((char*) &Eotvos, sizeof(double));
        file.read((char*) &resolution, sizeof(double));
        file.read((char*) &rho_l, sizeof(double));
        file.read((char*) &gamma, sizeof(double));
        // file.read((char*) &diameter, sizeof(double));
        file.read((char*) &mu_ratio, sizeof(double));
        // file.read((char*) &c_s, sizeof(double));
        // file.read((char*) &sigma, sizeof(double));
        // file.read((char*) &g, sizeof(double));
        file.read((char*) &s_3, sizeof(double));
        file.read((char*) &s_5, sizeof(double));

        file.read((char*) &width, sizeof(int));
        file.read((char*) &height, sizeof(int));

        Preprocess prepro(ReynoldsMax, Morton, Eotvos, resolution, rho_l, gamma, mu_ratio, s_3, s_5, width, height);
        
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

void write_techplot_output(const Lattice& l, int iterNum, bool verbose)
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
        PsiFile << "\"ux1\""<<endl;
        PsiFile << "\"uy1\""<<endl;
        PsiFile << "\"ux2\""<<endl;
        PsiFile << "\"uy2\""<<endl;
        PsiFile << "\"u_abs\""<<endl;
        PsiFile << "\"grad_x\""<<endl;
        PsiFile << "\"grad_y\""<<endl;
    }
    PsiFile << "ZONE T=\" "<< iterNum << "\""<< endl;
    PsiFile << "I="<<xsize<<", J="<<ysize<<", K=1,F=POINT"<< endl;
    if (verbose == true) PsiFile << "DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)"<<endl;
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
                VeloSet u = tmp.getU();
                ColSet rho = tmp.getRho();
                Vector v;
                Vector gradient = l.getGradient(i, j);
                if(sum(rho) > 0) v = (u[0]*rho[0] + u[1]*rho[1]) * (1/sum(rho));
                PsiFile << "\t" << sum(rho) << "\t" << u[0].x << "\t" << u[0].y << "\t" << u[1].x << "\t" << u[1].y << "\t" << v.Abs() << "\t" << gradient.x << "\t" << gradient.y;
            }
            PsiFile << endl;
        }
    }
    PsiFile.close();
}
void write_techplot_output_alternative(const Lattice& l, const string& filename)
{
    ofstream PsiFile;
    Cell tmp;
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
            
            VeloSet u = tmp.getU();
            ColSet rho = tmp.getRho();
            Vector v;
            Vector gradient = l.getGradient(i, j);
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

void write_vtk_output(const Lattice& l, int iterNum)
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

    VTKFile << "\nVECTORS u1 DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();

            VeloSet u = tmp.getU();
            VTKFile << u[0].x << " " << u[0].y << " 0 ";
        }
    }

    VTKFile << "\nVECTORS u2 DOUBLE"<<endl;
    for (int j = 0; j < ysize; j++)
    {
        for (int i = 0; i < xsize; i++)
        {
            tmp = l.getCell(i,j);
            tmp.calcRho();

            VeloSet u = tmp.getU();
            VTKFile << u[1].x << " " << u[1].y << " 0 ";
        }
    }
    VTKFile.close();
}

void write_data_plot(const std::vector<double> x, const std::vector<double> y1, const std::vector<double> y2, const string& filename){
    ofstream massplot;

    stringstream name;
    name << filename;

    massplot.open( name.str().c_str() );
    massplot << "x" << "\t" << "y1" << "\t" << "y2" << "\n";
    for(unsigned int i = 0; i< x.size(); i++){
        massplot << x.at(i) << "\t" << y1.at(i) << "\t" << y2.at(i) << "\n";
    } 
    massplot.close();
}

void write_data_plot(const std::vector<double> y, double del_x, const string& filename){
    ofstream RePlot;

    stringstream name;
    name << filename;

    RePlot.open( name.str().c_str() );
    RePlot << "x" << "\t" << "y" << "\n";
    for(unsigned int i = 0; i< y.size(); i++){
        RePlot << i*del_x << "\t" << y.at(i) << "\n";
    } 
    RePlot.close();
}

void write_param_log(const ParamSet& p){
    ofstream paramLog;

    stringstream name;
    name <<"paramLog";

    paramLog.open( name.str().c_str() );
    paramLog << "# used setup parameters" << endl;
    paramLog << "\n# omega [-]" << endl;
    paramLog << "omega_red = "  << p.getOmegaRed()           << endl;
    paramLog << "omega_blue = " << p.getOmegaBlue()          << endl;
    paramLog << "\n# density of denser phase, normalized [-]" << endl;
    paramLog << "rho_red = "    << p.getRhoR()               << endl;
    paramLog << "\n# density ratio [-]" << endl;
    paramLog << "gamma = "      << p.getGamma()              << endl;
    paramLog << "\n# method specific parameters [-]" << endl;
    paramLog << "alpha_blue = " << p.getAlpha()              << endl;
    paramLog << "delta = "      << p.getInterfaceThickness() << endl;
    paramLog << "beta = "       << p.getBeta()               << endl;
    paramLog << "\n# surface tension [-]" << endl;
    paramLog << "sigma = "      << p.getSigma()              << endl; //TODO get SI units of sigma
    paramLog << "\n# discretisization [-]" << endl;
    paramLog << "speedlimit = " << p.getSpeedlimit()         << endl; 
    paramLog << "dt = "         << p.getDeltaT()             << endl;
    paramLog << "dx = "         << p.getDeltaX()             << endl;
    paramLog << "gravity = "    << p.getG()                  << endl;

    RelaxationPar relax = p.getRelaxation();
    paramLog << "\n# MRT parameters [-]" << endl;
    paramLog << "s_2 = " << relax.s_2 << endl;
    paramLog << "s_3 = " << relax.s_3 << endl;
    paramLog << "s_5 = " << relax.s_5 << endl;


    paramLog.close();
}

//=========================== READ INPUT ===========================

const Preprocess read_preprocess_file(const string& filename){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
        mm.insert(pair<string,double>("Reynolds",10));
        mm.insert(pair<string,double>("Morton",100));
        mm.insert(pair<string,double>("Eotvos",10));
        mm.insert(pair<string,double>("resolution",30));
        mm.insert(pair<string,double>("rho_l",1));
        mm.insert(pair<string,double>("gamma",2));
        // mm.insert(pair<string,double>("diameter",0.1));
        // mm.insert(pair<string,double>("c_s",10));
        // mm.insert(pair<string,double>("sigma",1e-4));
        // mm.insert(pair<string,double>("g",10));
        mm.insert(pair<string,double>("mu_ratio",2));
        mm.insert(pair<string,double>("s_3",1));
        mm.insert(pair<string,double>("s_5",1));
        mm.insert(pair<string,double>("width",120));
        mm.insert(pair<string,double>("height",360));

        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
        }
    Preprocess prepro(mm.at("Reynolds"),mm.at("Morton"),mm.at("Eotvos"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"), mm.at("mu_ratio"), mm.at("s_3"), mm.at("s_5"), mm.at("width"), mm.at("height"));
    return prepro;
}

const Timetrack read_timetrack_file(const Preprocess& prepro, const string& filename){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("factor",1.1));
    mm.insert(pair<string,double>("max_steps",1e5));
    mm.insert(pair<string,double>("techplot_interval",1e3));
    mm.insert(pair<string,double>("restart_interval",1e4));
    
    // cycling through the input file
    double tmp;
    for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
        if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
    }

    Timetrack time(prepro.getTimestep(), mm.at("factor"), mm.at("max_steps"), mm.at("techplot_interval"), mm.at("restart_interval"));
 
    return time;
}

const ParamSet read_paramset_file(const string& filename){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
        mm.insert(pair<string,double>("omega_red",1));
        mm.insert(pair<string,double>("omega_blue",1));
        mm.insert(pair<string,double>("rho_red",1));
        mm.insert(pair<string,double>("gamma",2));
        mm.insert(pair<string,double>("sigma",1e-4));
        mm.insert(pair<string,double>("gravity",1e-4));
        mm.insert(pair<string,double>("speedlimit",1));
        mm.insert(pair<string,double>("dt",1e-3));
        mm.insert(pair<string,double>("dx",1e-3));

        mm.insert(pair<string,double>("s_2",1));
        mm.insert(pair<string,double>("s_3",1));
        mm.insert(pair<string,double>("s_5",1));

        mm.insert(pair<string,double>("alpha_blue",0.2));
        mm.insert(pair<string,double>("delta",0.1));
        mm.insert(pair<string,double>("beta",0.99));
     

        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
        }

    RelaxationPar rel(mm.at("s_2"),mm.at("s_3"),mm.at("s_5"));

    ParamSet params(mm.at("omega_red"), mm.at("omega_blue"), mm.at("rho_red"), mm.at("gamma"), mm.at("sigma"), mm.at("gravity"), mm.at("speedlimit"), mm.at("dt"), mm.at("dx"), rel, mm.at("alpha_blue"), mm.at("delta"), mm.at("beta")); /// < consructor
    return params;
}

const bool input_query(const string& filename, const string& query, double& value){
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

//=========================== AUXILIARY ===========================

const string createFilename(const string& name, int iteration, const string& type)
{
    stringstream filename;
    filename << name << iteration << type;
    return filename.str();
}