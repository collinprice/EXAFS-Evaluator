#include "dcdhelper.h"

extern "C" {
#include "dcd.h"
}

#include <exception>
#include <iostream>

#define C_TEXT( text ) ((char*)std::string( text ).c_str())

std::vector< std::vector<PDBAtom> > DCDHelper::getXYZs(std::string filename, int size) {

	initialize_read_xyz(C_TEXT(filename));
	
	float* x = new float[my_H.N];
	float* y = new float[my_H.N];
	float* z = new float[my_H.N];

	if (size < 1) {
		throw "Size must be greater than zero.";
	} else if (size > my_H.NSet) {
		throw "DCD File does not have enough frames for requested size.";
	}

	double percent = (double)size / my_H.NSet;

	int step = my_H.NSet / (my_H.NSet * percent);
	std::vector< std::vector<PDBAtom> > molecules;
	int counter = 0;
	for (int j = 0; j < my_H.NSet; j += step) {
		++counter;

		read_xyz(j, x, y, z);
		std::vector<PDBAtom> points;
		for (int i = 0; i < my_H.N; ++i) {
			points.push_back(PDBAtom(i, x[i], y[i], z[i]));
		}
		molecules.push_back(points);
	}
	
	delete[] x;
	delete[] y;
	delete[] z;

	return molecules;
}