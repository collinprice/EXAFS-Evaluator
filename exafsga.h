#include "exafsevaluator.h"
#include "chromosome.h"

#include <vector>
#include <fstream>

class EXAFSGA {
	public:

		EXAFSGA(EXAFSEvaluator* exafs_evaluator, double mutation_rate, double crossover_rate, bool elitism, int max_generations, std::string results_file);
		~EXAFSGA();
		void begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations);
		void begin_recentering(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, int population_size, double convergence_rate, int max_iterations);

	private:

		std::vector<Chromosome> population;
		EXAFSEvaluator* exafs_evaluator;
		double mutation_rate;
		double crossover_rate;
		bool elitism;
		int max_generations;
		std::string results_file;
		std::ofstream output_stream;
		std::vector<Chromosome> best_individuals;
		std::string stats_folder;

		void initPopulation(std::vector< std::vector<PDBAtom> > population);

		void evolve();
		Chromosome selection();

		void evaluate( Chromosome& child );
		void evaluatePopulation();

		void crossover(Chromosome& p1, Chromosome& p2);
		void mutate(Chromosome& child);

		Chromosome best_chromosome();
		void saveBestChromosome();

		double unifRand();

		void initStats();
		void recordStats();
		void finalStats();

		bool convergence();
		bool convergence(double rate);
};