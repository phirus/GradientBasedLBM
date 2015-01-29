#include"BasicIO.h"


//=========================== WRITE OUTPUT ===========================

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

    RePlot.open( filename.c_str() );
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
    paramLog << "dt = "         << p.getDeltaT()             << endl;
    paramLog << "dx = "         << p.getDeltaX()             << endl;
    paramLog << "gravity = "    << p.getG()                  << endl;

    RelaxationPar3D relax = p.getRelaxation3D();
    paramLog << "\n# MRT parameters [-]" << endl;
    paramLog << "s_2 = " << relax.s_2 << endl;
    paramLog << "s_3 = " << relax.s_3 << endl;
    paramLog << "s_5 = " << relax.s_5 << endl;
    paramLog << "s_11 = " << relax.s_11 << endl;
    paramLog << "s_17 = " << relax.s_17 << endl;


    paramLog.close();
}


//=========================== READ INPUT ===========================

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
        mm.insert(pair<string,double>("mu_ratio",2));
        mm.insert(pair<string,double>("s_3",1));
        mm.insert(pair<string,double>("s_5",1));
        mm.insert(pair<string,double>("s_11",1));
        mm.insert(pair<string,double>("s_17",1));
        mm.insert(pair<string,double>("xCells",50));
        mm.insert(pair<string,double>("yCells",50));
        mm.insert(pair<string,double>("zCells",50));

        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
        }
    Preprocess prepro(mm.at("Reynolds"),mm.at("Morton"),mm.at("Eotvos"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"), mm.at("mu_ratio"), mm.at("s_3"), mm.at("s_5"), mm.at("s_11"), mm.at("s_17"), mm.at("xCells"), mm.at("yCells"), mm.at("zCells"));
    return prepro;
}

const Timetrack read_timetrack_file(const string& filename){
    vector<string> tags;
    vector<double> val;
    // initialzing strings and fallback values
    map<string,double> mm;
    mm.insert(pair<string,double>("max_steps",1e5));
    mm.insert(pair<string,double>("output_interval",2e3));
    mm.insert(pair<string,double>("restart_interval",1e4));
    
    // cycling through the input file
    double tmp;
    for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
        if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
    }

    Timetrack time(mm.at("max_steps"), mm.at("output_interval"), mm.at("restart_interval"));
 
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
     

        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
        }

    const RelaxationPar3D rel(mm.at("s_2"),mm.at("s_3"),mm.at("s_5"),mm.at("s_11"),mm.at("s_17"));

    ParamSet params(mm.at("omega_red"), mm.at("omega_blue"), mm.at("rho_red"), mm.at("gamma"), mm.at("sigma"), mm.at("gravity"), mm.at("dt"), mm.at("dx"), rel, mm.at("alpha_blue"), mm.at("delta"), mm.at("beta")); /// < consructor
    return params;
}


//=========================== AUXILIARY ===========================

const string createFilename(const string& name, int iteration, const string& type)
{
    stringstream filename;
    filename << name << iteration << type;
    return filename.str();
}