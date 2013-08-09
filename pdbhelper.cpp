#include "pdbhelper.h"

extern "C" {
#include "file_read_write.h"
}

#include <iostream>
#include <fstream>

#define C_TEXT( text ) ((char*)std::string( text ).c_str())

PDBHelper::PDBHelper(std::string pdb_file, std::string amber_topology_file, std::string output_pdb_file) {

	read_pdb(C_TEXT(pdb_file));
	read_parm7(C_TEXT(amber_topology_file)); // Restore reliable atom names
	mass_to_element();

	this->output_pdb_file = output_pdb_file;
}

PDBHelper::~PDBHelper() {

	std::ifstream temp_file(this->output_pdb_file.c_str());
	if (temp_file.good()) {
		system(("rm " + this->output_pdb_file).c_str());
	}
}

int PDBHelper::numberOfAtoms() {
	return N;
}

std::string PDBHelper::atomAtIndex(int index) {
	return std::string(1, Element[index][0]) + std::string(1, Element[index][1]);
}

std::vector<PDBAtom> PDBHelper::getEXAFSAtoms() {

	std::vector<PDBAtom> exafs_atoms;

	for (int i = 0; i < N; ++i) {
		if (occupancy[i] == 1) {
			exafs_atoms.push_back(PDBAtom(this->atomAtIndex(i), i, X[i], Y[i], Z[i]));
		}
	}

	return exafs_atoms;
}
	
void PDBHelper::updateEXAFSAtoms(std::vector<PDBAtom> atoms) {

	for (std::vector<PDBAtom>::iterator atom = atoms.begin(); atom != atoms.end(); ++atom) {
		
		int index = atom->getIndex();
		X[index] = atom->x;
		Y[index] = atom->y;
		Z[index] = atom->z;
	}
}

void PDBHelper::updateAtomsFromXYZ(std::string filename, bool onlyEXAFS) {

	std::ifstream xyz_file(filename.c_str());
	std::string x,y,z;

	int index = -1;
	while(xyz_file.good()) {

		xyz_file >> x >> y >> z;
		++index;

		if (onlyEXAFS && occupancy[index] != 1) continue;

		X[index] = atof(x.c_str());
		Y[index] = atof(y.c_str());
		Z[index] = atof(z.c_str());		
	}

	xyz_file.close();
}

void PDBHelper::writePDBFile(std::string filename) {
	write_one_pdb(C_TEXT(filename));
}

void PDBHelper::writePDBFile() {
	this->writePDBFile(this->output_pdb_file);
}