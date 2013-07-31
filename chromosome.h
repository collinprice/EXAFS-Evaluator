#include "pdbatom.h"

#include <vector>

class Chromosome: public std::vector<PDBAtom> {

	public:
		double exafs_score;
		double potential_energy;

		Chromosome();
		Chromosome(std::vector<PDBAtom> points);
};