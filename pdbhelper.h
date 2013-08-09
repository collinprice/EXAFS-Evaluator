#include "pdbatom.h"

#include <string>
#include <vector>

class PDBHelper {

public:

	int numberOfAtoms();
	std::string atomAtIndex(int index);
	std::string output_pdb_file;

	PDBHelper(std::string pdb_file, std::string amber_topology_file, std::string output_pdb_file);
	~PDBHelper();

	std::vector<PDBAtom> getEXAFSAtoms();
	void updateEXAFSAtoms(std::vector<PDBAtom> atoms);
	void updateAtomsFromXYZ(std::string filename, bool onlyEXAFS);
	void writePDBFile(std::string filename);
	void writePDBFile();

};