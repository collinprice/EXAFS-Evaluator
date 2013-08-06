#include "pdbatom.h"

#include <vector>

class Chromosome {

	public:
		double exafs_score;
		double potential_energy;
		bool is_evaluated;

		std::vector< std::pair<double, double> > exafs_data;
		std::vector< PDBAtom > atoms;

		Chromosome();
		Chromosome( std::vector<PDBAtom> atoms );
		Chromosome( const Chromosome& other );

	private:

		void init();
};