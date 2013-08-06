#include "exafsevaluator.h"
#include "chromosome.h"

#include <vector>
#include <fstream>

class EXAFSGA {
	public:

		EXAFSGA(EXAFSEvaluator* exafs_evaluator, double mutation_rate, double crossover_rate, bool elitism, int max_generations, std::string results_file);
		~EXAFSGA();
		void begin(std::vector< std::vector<PDBAtom> > initial_population);

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

		void initPopulation(std::vector< std::vector<PDBAtom> > population);

		void evolve();
		Chromosome selection();

		void evaluate( Chromosome& child );
		void evaluatePopulation();

		void crossover(Chromosome& p1, Chromosome& p2);
		void mutate(Chromosome& child);

		Chromosome best_chromosome();

		double unifRand();

		void initStats();
		void recordStats();
		void finalStats();
};