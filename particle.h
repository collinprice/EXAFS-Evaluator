#include "chromosome.h"

class Particle : public Chromosome {

	public:
		double best_exafs_score;
		double best_potential_energy;

		std::vector< std::pair<double, double> > best_exafs_data;
		std::vector< PDBAtom > best_atoms;

		std::vector< PDBAtom > velocity;

		Particle();
		Particle(double range, std::vector<PDBAtom> atoms);
		Particle(const Particle& other);

		void updateBest();
		void updatePosition();
		void updateVelocity(Particle global_best, double inertia, double social, double cognitive);
	
	private:
		double unifRand();
		double randomRange(double max);
};