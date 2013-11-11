#include "exafsevaluator.h"
#include "particle.h"

#include <vector>
#include <fstream>

class EXAFSPSO {
	public:

		EXAFSPSO(EXAFSEvaluator* exafs_evaluator, double mutation_rate, double crossover_rate, bool elitism, int max_generations, std::string results_file);
		~EXAFSPSO();
		void begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, double max_initial_velocty_offset);

	private:

		std::vector<Particle> population;
		EXAFSEvaluator* exafs_evaluator;
		double mutation_rate;
		double crossover_rate;
		bool elitism;
		int max_generations;
		std::string results_file;
		std::ofstream output_stream;
		std::vector<Particle> best_individuals;
		std::string stats_folder;
		Particle best_particle;

		void initPopulation(std::vector< std::vector<PDBAtom> > population, double max_initial_velocty_offset);

		void evolve();

		void evaluate( Particle& child );
		void evaluatePopulation();

		Particle best_chromosome();
		void saveBestChromosome();

		double unifRand();

		void initStats();
		void recordStats();
		void finalStats();

		bool convergence();
		bool convergence(double rate);
};