#include "exafsevaluator.h"
#include "chromosome.h"

#include <vector>
#include <fstream>

class EXAFSDE {
	public:

		EXAFSDE(EXAFSEvaluator* exafs_evaluator, double f, double cr, int max_generations, std::string results_file);
		~EXAFSDE();
		void begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations);
		void begin_recentering(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, int population_size, double keep_percentage, int max_iterations);

	private:

		std::vector<Chromosome> population;
		EXAFSEvaluator* exafs_evaluator;
		double f;
		double cr;
		int max_generations;
		std::string results_file;
		std::ofstream output_stream;
		std::vector<Chromosome> best_individuals;
		std::string stats_folder;

		void initPopulation(std::vector< std::vector<PDBAtom> > population);

		void evolve();

		Chromosome mutate(int i, int r1, int r2, int r3);

		void evaluate( Chromosome& child );
		void evaluateRMSD( Chromosome& child );
		void evaluatePotentialEnergy( Chromosome& child );
		
		void evaluatePopulation();

		Chromosome best_chromosome();
		void saveBestChromosome();

		double unifRand();
		double randInt(int max);

		void initStats();
		void recordStats();
		void finalStats();

		bool convergence();
		bool convergence(double rate);
};