#include "pdbatom.h"

#include <string>
#include <vector>

class DCDHelper {

public:

	static std::vector< std::vector<PDBAtom> > getXYZs(std::string filename);

};