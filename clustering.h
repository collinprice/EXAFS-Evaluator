#include "pdbatom.h"

#include <vector>

class Clustering {

public:

	static std::vector< std::vector<double> > createTable(std::vector< std::vector<PDBAtom> > atoms);

};