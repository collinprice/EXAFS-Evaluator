#include "chromosome.h"

class Particle : public Chromosome {

	public:
		double best_exafs_score;
		double best_potential_energy;

		std::vector< std::pair<double, double> > best_exafs_data;
		std::vector< PDBAtom > best_atoms;

		std::vector< PDBAtom > velocity;

		Particle();
	
};