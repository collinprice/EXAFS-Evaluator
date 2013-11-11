#include "pdbatom.h"

#include <vector>

class Particle {

	public:

		double exafs_score;
		double potential_energy;
		bool is_evaluated;

		std::vector< std::pair<double, double> > exafs_data;
		std::vector< PDBAtom > position;
		std::vector< std::vector<double> > velocity;

		std::vector< PDBAtom > best_position;
		double best_exafs_score;

		Particle();
		Particle( std::vector<PDBAtom> atoms, double max_offset );
		Particle( const Particle& other );

		void update_velocity(std::vector< PDBAtom > global_best_position);
		void update_position();
		void update_best_position();

		bool operator==(const Particle& c) { return this->exafs_score == c.exafs_score; };

	private:

		void init();
		void generate_initial_velocity(double max_offset);
};