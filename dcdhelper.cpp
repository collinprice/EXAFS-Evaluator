#include "dcdhelper.h"

extern "C" {
#include "dcd.h"
}

#include <iostream>

#define C_TEXT( text ) ((char*)std::string( text ).c_str())

std::vector< std::vector<PDBAtom> > DCDHelper::getXYZs(std::string filename) {

	initialize_read_xyz(C_TEXT(filename));
	
	float* x = new float[my_H.N];
	float* y = new float[my_H.N];
	float* z = new float[my_H.N];

	std::cout << "Number of frames: " << my_H.NSet << std::endl;

	std::vector< std::vector<PDBAtom> > molecules;
	for (int j = 0; j < my_H.NSet; j += 100) {
		
		std::cout << j << std::endl;

		read_xyz(j, x, y, z);
		std::vector<PDBAtom> points;
		for (int i = 0; i < my_H.N; ++i) {
			points.push_back(PDBAtom(i, x[i], y[i], z[i]));
		}
		molecules.push_back(points);
	}

	free(x);
	free(y);
	free(z);

	return molecules;
}