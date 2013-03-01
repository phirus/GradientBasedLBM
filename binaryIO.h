#ifndef BINARYIO_H
#define BINARYIO_H

#include<fstream>
#include<iostream>
#include"lattice.h"

void binary_output(const Lattice& l);
const Lattice binary_input();

#endif

