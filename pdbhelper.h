#include "pdbatom.h"

#include <string>
#include <vector>

class PDBHelper {

public:
	PDBHelper(std::string pdb_file, std::string amber_topology_file);

	int numberOfAtoms();
	std::string atomAtIndex(int index);

	std::vector<PDBAtom> getEXAFSAtoms();
	void updateEXAFSAtoms(std::vector<PDBAtom> atoms);
	void writePDBFile(std::string filename);
};