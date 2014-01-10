#include "exafsevaluator.h"
#include "particle.h"

#include <vector>
#include <fstream>

class EXAFSPSO {
	public:

		EXAFSPSO(EXAFSEvaluator* exafs_evaluator, double inertia, double social, double cognitive, double velocity_range, int max_generations, std::string results_file);
		~EXAFSPSO();
		void begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations);
		// void begin_recentering(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, int population_size, double keep_percentage, int max_iterations);

	private:

		std::vector<Particle> population;
		Particle global_best_particle;

		EXAFSEvaluator* exafs_evaluator;
		double inertia;
		double social;
		double cognitive;
		double velocity_range;
		
		int max_generations;
		std::string results_file;
		std::ofstream output_stream;
		std::vector<Particle> best_individuals;
		std::string stats_folder;

		void initPopulation(std::vector< std::vector<PDBAtom> > population);

		void updateVelocities();
		void updatePositions();
		void updateGlobalBest();

		void evaluate( Particle& child );
		void evaluatePopulation();

		Particle best_particle();
		void saveBestParticle();

		double unifRand();
		double randInt(int max);

		void initStats();
		void recordStats();
		void finalStats();

		bool convergence();
		bool convergence(double rate);
};