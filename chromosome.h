#pragma once
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

		bool operator==(const Chromosome& c) { return this->exafs_score == c.exafs_score; };

	private:

		void init();
};