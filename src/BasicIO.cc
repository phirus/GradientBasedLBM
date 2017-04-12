#include"BasicIO.h"


//=========================== WRITE OUTPUT ===========================

void write_file_header(const string& filename, const string& header)
{
    ofstream Plot;
    Plot.open( filename.c_str() );
    Plot << header << "\n";
    Plot.close();
}

void write_data_plot(const std::vector<double> y, double del_x, const string& filename)
{
    ofstream RePlot;

    RePlot.open( filename.c_str() );
    RePlot << "x" << "\t" << "y" << "\n";
    for(unsigned int i = 0; i< y.size(); i++){
        RePlot << i*del_x << "\t" << y.at(i) << "\n";
    } 
    RePlot.close();
}

void write_data_plot_linewise(int time ,double y1, double y2, const string& filename)
{
    ofstream Plot;
    Plot.open( filename.c_str(), ios::app);
    Plot << time << "\t" << y1 << "\t" << y2 << "\n"; 
    Plot.close();
}

void write_csv(const nested_vector& data, const string& filename, const string& header)
{
    const auto numOfSets = data.size();
    const auto numOfLines = data.at(0).size();

    ofstream Plot;
    Plot.open( filename.c_str() );
    Plot << header << "\n";

    for (unsigned int line = 0; line  < numOfLines ; line++)
    {
        for (unsigned int item = 0; item < numOfSets; item++)
        {
            const string end_char = (item == numOfSets -1) ? "\n" : ";" ;
            Plot << data[item][line] << end_char;
        }
    } 
    Plot.close();
}

void write_csv_linewise(int i, double Posx, double PosY, double v_x, double v_y, double Re, double F_x, double F_y, const string& filename)
{
    ofstream Plot;
    Plot.open(filename.c_str(),ios::app);
    Plot << i << ";" << Posx << ";" << PosY << ";" << v_x << ";" << v_y << ";" << Re << ";" << F_x << ";" << F_y << "\n";
    Plot.close();
}

void write_param_log(const ParamSet& p)
{
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
    paramLog << "dt = "         << p.getDeltaT()             << endl;
    paramLog << "dx = "         << p.getDeltaX()             << endl;
    paramLog << "gravity = "    << p.getG()                  << endl;
    paramLog << "bulk_visco = " << p.getBulkVisco()          << endl;

    RelaxationPar3D relax = p.getRelaxation3D();
    paramLog << "\n# MRT parameters [-]" << endl;
    paramLog << "s_2 = " << relax.s_2 << endl;
    paramLog << "s_3 = " << relax.s_3 << endl;
    paramLog << "s_5 = " << relax.s_5 << endl;
    paramLog << "s_11 = " << relax.s_11 << endl;
    paramLog << "s_17 = " << relax.s_17 << endl;

    paramLog.close();
}

void write_param_log_csv(const ParamSet& p)
{
    ofstream paramLog;

    stringstream name;
    name <<"paramLog.csv";
    //paramLog.precision(std::numeric_limits< double >::max_digits10);
    paramLog.open( name.str().c_str() );

    paramLog << "omega_red;omega_blue;rho_red;gamma;alpha_blue;delta;beta;sigma;dt;dx;gravity;bulk_visco;s_2;s_3;s_5;s_11;s_17" << endl;
    RelaxationPar3D relax = p.getRelaxation3D();
    paramLog <<  p.getOmegaRed() << ";" << p.getOmegaBlue()  << ";" << p.getRhoR() << ";" << p.getGamma() << ";" ;
    paramLog << p.getAlpha() << ";" << p.getInterfaceThickness() << ";" << p.getBeta() << ";" << p.getSigma() << ";" ; 
    paramLog << p.getDeltaT() << ";" << p.getDeltaX() << ";" << p.getG() << ";" << p.getBulkVisco() << ";";
    paramLog << relax.s_2 << ";" << relax.s_3 << ";" << relax.s_5 << ";" << relax.s_11 << ";" << relax.s_17 << endl;

    paramLog.close();
}

void write_preprocess_csv(const Preprocess& p)
{
    ofstream preproCSV;

    stringstream name;
    name <<"prepro.csv";
    preproCSV.open( name.str().c_str() );

    preproCSV << "ReynoldsMax;Morton;Eotvos;resolution;muRatio;bulk_visco;isShearFlow;shearRate;xCells;yCells;zCells"<<endl;
    preproCSV << p.getReynoldsMax() << ";" << p.getMorton() << ";" << p.getEotvos() << ";" << p.getResolution() << ";" ;
    preproCSV << p.getMuRatio() << ";" << p.getBulkVisco() << ";" << p.getIsShearFlow() << ";" << p.getShearRate() << ";" ;
    preproCSV << p.getXCells() << ";" << p.getYCells() << ";" << p.getZCells() << endl;
    preproCSV.close();
}

void write_preprocess_drop_csv(const Preprocess_Drop& p)
{
    ofstream preproCSV;

    stringstream name;
    name <<"prepro.csv";
    preproCSV.open( name.str().c_str() );

    preproCSV << "VeloFrac;Ohnesorge;Weber;resolution;muRatio;bulk_visco;xCells;yCells;zCells"<<endl;
    preproCSV << p.getVeloFrac() << ";" << p.getOhnesorge() << ";" << p.getWeber() << ";" << p.getResolution() << ";" ;
    preproCSV << p.getMuRatio() << ";" << p.getBulkVisco() << ";" ;
    preproCSV << p.getXCells() << ";" << p.getYCells() << ";" << p.getZCells() << endl;
    preproCSV.close();
}

//=========================== READ INPUT ===========================

const bool input_query(const string& filename, const string& query, double& value)
{
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

const map<string,double> assign_map_via_file(map<string,double> mm, const string& filename)
{
    double tmp;
    for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){
        if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
    }
    return mm;
}

const Preprocess read_preprocess_file(const string& filename)
{

    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("Reynolds",10));
    mm.insert(pair<string,double>("Morton",100));
    mm.insert(pair<string,double>("Eotvos",10));
    mm.insert(pair<string,double>("resolution",30));
    mm.insert(pair<string,double>("rho_l",1));
    mm.insert(pair<string,double>("gamma",2));
    mm.insert(pair<string,double>("mu_ratio",1));
    mm.insert(pair<string,double>("bulk_visco",2));
    mm.insert(pair<string,double>("s_3",1));
    mm.insert(pair<string,double>("s_5",1));
    mm.insert(pair<string,double>("s_11",1));
    mm.insert(pair<string,double>("s_17",1));
    mm.insert(pair<string,double>("isShearFlow",0)); // implicit type cast: 0 -> false, 1 -> true
    mm.insert(pair<string,double>("shearRate",0));
    mm.insert(pair<string,double>("xCells",50));
    mm.insert(pair<string,double>("yCells",50));
    mm.insert(pair<string,double>("zCells",50));

    mm = assign_map_via_file(mm, filename);

    Preprocess prepro(mm.at("Reynolds"),mm.at("Morton"),mm.at("Eotvos"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"), mm.at("mu_ratio"), mm.at("bulk_visco"), mm.at("s_3"), mm.at("s_5"), mm.at("s_11"), mm.at("s_17"), mm.at("isShearFlow"), mm.at("shearRate"), mm.at("xCells"), mm.at("yCells"), mm.at("zCells"));
    return prepro;
}

const Preprocess_Drop read_preprocess_drop_file(const string& filename)
{

    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("VeloFrac",10));
    mm.insert(pair<string,double>("Ohnesorge",100));
    mm.insert(pair<string,double>("Weber",10));
    mm.insert(pair<string,double>("resolution",30));
    mm.insert(pair<string,double>("rho_l",1));
    mm.insert(pair<string,double>("gamma",2));
    mm.insert(pair<string,double>("mu_ratio",1));
    mm.insert(pair<string,double>("bulk_visco",2));
    mm.insert(pair<string,double>("s_3",1));
    mm.insert(pair<string,double>("s_5",1));
    mm.insert(pair<string,double>("s_11",1));
    mm.insert(pair<string,double>("s_17",1));
    mm.insert(pair<string,double>("xCells",50));
    mm.insert(pair<string,double>("yCells",50));
    mm.insert(pair<string,double>("zCells",50));

    mm = assign_map_via_file(mm, filename);

    Preprocess_Drop prepro(mm.at("VeloFrac"),mm.at("Ohnesorge"),mm.at("Weber"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"), mm.at("mu_ratio"), mm.at("bulk_visco"), mm.at("s_3"), mm.at("s_5"), mm.at("s_11"), mm.at("s_17"), mm.at("xCells"), mm.at("yCells"), mm.at("zCells"));
    return prepro;
}

const Timetrack read_timetrack_file(const string& filename)
{
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("max_steps",1e5));
    mm.insert(pair<string,double>("output_interval",2e3));
    mm.insert(pair<string,double>("restart_interval",1e4));
    
    mm = assign_map_via_file(mm, filename);

    Timetrack time(mm.at("max_steps"), mm.at("output_interval"), mm.at("restart_interval"));
 
    return time;
}

const ParamSet read_paramset_file(const string& filename)
{
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("omega_red",1));
    mm.insert(pair<string,double>("omega_blue",1));
    mm.insert(pair<string,double>("rho_red",1));
    mm.insert(pair<string,double>("gamma",2));
    mm.insert(pair<string,double>("sigma",1e-4));
    mm.insert(pair<string,double>("gravity",1e-4));
    mm.insert(pair<string,double>("bulk_visco",2));
    mm.insert(pair<string,double>("dt",1e-3));
    mm.insert(pair<string,double>("dx",1e-3));
    mm.insert(pair<string,double>("s_2",1));
    mm.insert(pair<string,double>("s_3",1));
    mm.insert(pair<string,double>("s_5",1));
    mm.insert(pair<string,double>("s_11",1));
    mm.insert(pair<string,double>("s_17",1));
    mm.insert(pair<string,double>("alpha_blue",0.2));
    mm.insert(pair<string,double>("delta",0.1));
    mm.insert(pair<string,double>("beta",0.99));
     
    mm = assign_map_via_file(mm, filename);

    const RelaxationPar3D rel(mm.at("omega_red"),mm.at("s_2"),mm.at("s_3"),mm.at("s_5"),mm.at("s_11"),mm.at("s_17"));

    ParamSet params(mm.at("omega_red"), mm.at("omega_blue"), mm.at("rho_red"), mm.at("gamma"), mm.at("sigma"), mm.at("gravity"), mm.at("dt"), mm.at("dx"), rel, mm.at("bulk_visco"), mm.at("alpha_blue"), mm.at("delta"), mm.at("beta")); /// < consructor
    return params;
}

const Boundaries read_boundaries_file(const string& filename)
{   
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("n_type",0));
    mm.insert(pair<string,double>("n_rho_dense",0));
    mm.insert(pair<string,double>("n_rho_dilute",0));
    mm.insert(pair<string,double>("nvx",0));
    mm.insert(pair<string,double>("nvy",0));
    mm.insert(pair<string,double>("nvz",0));

    mm.insert(pair<string,double>("s_type",0));
    mm.insert(pair<string,double>("s_rho_dense",0));
    mm.insert(pair<string,double>("s_rho_dilute",0));
    mm.insert(pair<string,double>("svx",0));
    mm.insert(pair<string,double>("svy",0));
    mm.insert(pair<string,double>("svz",0));

    mm.insert(pair<string,double>("e_type",0));
    mm.insert(pair<string,double>("e_rho_dense",0));
    mm.insert(pair<string,double>("e_rho_dilute",0));
    mm.insert(pair<string,double>("evx",0));
    mm.insert(pair<string,double>("evy",0));
    mm.insert(pair<string,double>("evz",0));

    mm.insert(pair<string,double>("w_type",0));
    mm.insert(pair<string,double>("w_rho_dense",0));
    mm.insert(pair<string,double>("w_rho_dilute",0));
    mm.insert(pair<string,double>("wvx",0));
    mm.insert(pair<string,double>("wvy",0));
    mm.insert(pair<string,double>("wvz",0));

    mm.insert(pair<string,double>("f_type",0));
    mm.insert(pair<string,double>("f_rho_dense",0));
    mm.insert(pair<string,double>("f_rho_dilute",0));
    mm.insert(pair<string,double>("fvx",0));
    mm.insert(pair<string,double>("fvy",0));
    mm.insert(pair<string,double>("fvz",0));

    mm.insert(pair<string,double>("b_type",0));
    mm.insert(pair<string,double>("b_rho_dense",0));
    mm.insert(pair<string,double>("b_rho_dilute",0));
    mm.insert(pair<string,double>("bvx",0));
    mm.insert(pair<string,double>("bvy",0));
    mm.insert(pair<string,double>("bvz",0));

    mm = assign_map_via_file(mm, filename);

    ColSet tmp_rho;
    VeloSet3D tmp_u;

    Boundaries b;

    tmp_rho = {{mm.at("n_rho_dense"),mm.at("n_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("nvx"),mm.at("nvy"),mm.at("nvz")),Vector3D()}};
    b.north = BoundaryInformation(mm.at("n_type"),tmp_rho,tmp_u);

    tmp_rho = {{mm.at("s_rho_dense"),mm.at("s_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("svx"),mm.at("svy"),mm.at("svz")),Vector3D()}};
    b.south = BoundaryInformation(mm.at("s_type"),tmp_rho,tmp_u);

    tmp_rho = {{mm.at("e_rho_dense"),mm.at("e_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("evx"),mm.at("evy"),mm.at("evz")),Vector3D()}};
    b.east = BoundaryInformation(mm.at("e_type"),tmp_rho,tmp_u);

    tmp_rho = {{mm.at("w_rho_dense"),mm.at("w_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("wvx"),mm.at("wvy"),mm.at("wvz")),Vector3D()}};
    b.west = BoundaryInformation(mm.at("w_type"),tmp_rho,tmp_u);

    tmp_rho = {{mm.at("f_rho_dense"),mm.at("f_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("fvx"),mm.at("fvy"),mm.at("fvz")),Vector3D()}};
    b.front = BoundaryInformation(mm.at("f_type"),tmp_rho,tmp_u);

    tmp_rho = {{mm.at("b_rho_dense"),mm.at("b_rho_dilute")}};
    tmp_u = {{Vector3D(mm.at("bvx"),mm.at("bvy"),mm.at("bvz")),Vector3D()}};
    b.back = BoundaryInformation(mm.at("b_type"),tmp_rho,tmp_u);
 
    return b;
}

//=========================== AUXILIARY ===========================

const string createFilename(const string& name, int iteration, const string& type)
{
    stringstream filename;
    filename << name << iteration << type;
    return filename.str();
}