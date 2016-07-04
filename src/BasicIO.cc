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

void write_data_plot(const std::vector<double> y1, const std::vector<double> y2, double del_x, const string& filename){
    ofstream Plot;

    Plot.open( filename.c_str() );
    Plot << "x" << "\t" << "y1" << "\t" << "y2" << "\n";
    for(unsigned int i = 0; i< y1.size(); i++){
        Plot << i*del_x << "\t" << y1.at(i) << "\t" << y2.at(i) << "\n";
    } 
    Plot.close();
}

void write_csv(const nested_vector& data, const string& filename, const string& header){
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

void write_param_log_csv(const ParamSet& p){
    ofstream paramLog;

    stringstream name;
    name <<"paramLog.csv";
    //paramLog.precision(std::numeric_limits< double >::max_digits10);
    paramLog.open( name.str().c_str() );

    paramLog << "omega_red;omega_blue;rho_red;gamma;alpha_blue;delta;beta;sigma;dt;dx;gravity;s_2;s_3;s_5;s_11;s_17" << endl;
    RelaxationPar3D relax = p.getRelaxation3D();
    paramLog <<  p.getOmegaRed() << ";" << p.getOmegaBlue()  << ";" << p.getRhoR() << ";" << p.getGamma() << ";" ;
    paramLog << p.getAlpha() << ";" << p.getInterfaceThickness() << ";" << p.getBeta() << ";" << p.getSigma() << ";" ; 
    paramLog << p.getDeltaT() << ";" << p.getDeltaX() << ";" << p.getG() << ";";
    paramLog << relax.s_2 << ";" << relax.s_3 << ";" << relax.s_5 << ";" << relax.s_11 << ";" << relax.s_17 << endl;

    paramLog.close();
}

void write_preprocess_csv(const Preprocess& p){
    ofstream preproCSV;

    stringstream name;
    name <<"prepro.csv";
    preproCSV.open( name.str().c_str() );

    preproCSV << "ReynoldsMax;Morton;Eotvos;resolution;muRatio;isShearFlow;shearRate;xCells;yCells;zCells"<<endl;
    preproCSV << p.getReynoldsMax() << ";" << p.getMorton() << ";" << p.getEotvos() << ";" << p.getResolution() << ";" ;
    preproCSV << p.getMuRatio() << ";" << p.getIsShearFlow() << ";" << p.getShearRate() << ";" ;
    preproCSV << p.getXCells() << ";" << p.getYCells() << ";" << p.getZCells() << endl;
    preproCSV.close();
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
        mm.insert(pair<string,double>("isShearFlow",0)); // implicit type cast: 0 -> false, 1 -> true
        mm.insert(pair<string,double>("shearRate",0));
        mm.insert(pair<string,double>("xCells",50));
        mm.insert(pair<string,double>("yCells",50));
        mm.insert(pair<string,double>("zCells",50));

        // cycling through the input file
        double tmp;
        for(map<string,double>::iterator it = mm.begin(); it != mm.end(); it++){            
            if( input_query(filename,it->first,tmp) == true ) it->second = tmp;
        }
    Preprocess prepro(mm.at("Reynolds"),mm.at("Morton"),mm.at("Eotvos"),mm.at("resolution"),mm.at("rho_l"),mm.at("gamma"), mm.at("mu_ratio"), mm.at("s_3"), mm.at("s_5"), mm.at("s_11"), mm.at("s_17"), mm.at("isShearFlow"), mm.at("shearRate"), mm.at("xCells"), mm.at("yCells"), mm.at("zCells"));
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