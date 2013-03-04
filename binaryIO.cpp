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

    // setting up the file name
    stringstream name;
    name << filename;

    // setting up file
    fstream file(name.str().c_str(),ios::in | ios::binary);
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
