#include "pdbatom.h"

#include <string>
#include <vector>

class DCDHelper {

public:

	static std::vector< std::vector<PDBAtom> > getXYZs(std::string filename, int size);
	static std::vector< std::vector<PDBAtom> > getXYZs(std::string filename);
	static std::vector< std::vector<PDBAtom> > getXYZsByIndex(std::string filename, std::vector<int> indexes);

};