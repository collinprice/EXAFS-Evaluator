#include "exafsde.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

EXAFSDE::EXAFSDE(EXAFSEvaluator* exafs_evaluator, double f, double cr, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->f = f;
	this->cr = cr;
	this->max_generations = max_generations;
	this->results_file = results_file;
}

EXAFSDE::~EXAFSDE() {

	delete this->exafs_evaluator;
}

void EXAFSDE::begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {
		
		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		this->initPopulation(initial_populations[i]);
		this->initStats();
		this->recordStats();

		std::cout << "Begin Run " << (i+1) << std::endl;
		for (int i = 0; i < this->max_generations; ++i) {
			std::cout << "Generation: " << (i+1) << std::endl;
			this->evolve();

			this->saveBestChromosome();

			if (this->convergence()) break;
		}

		this->finalStats();
	}
}

bool chromosome_sort_de(Chromosome const & a, Chromosome const & b) {
	return a.exafs_score < b.exafs_score;
}

void EXAFSDE::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	this->best_individuals.clear();
	this->population.clear();
	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Chromosome child(*i);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
	this->saveBestChromosome();
}

void EXAFSDE::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		if (!this->population[i].is_evaluated) {
			this->evaluate(this->population[i]);
			this->population[i].is_evaluated = true;
		}
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << std::endl;
	}
}

void EXAFSDE::evaluate( Chromosome& child ) {

	this->exafs_evaluator->updateAtoms(child.atoms);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	// child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

void EXAFSDE::evolve() {
	
	std::vector<Chromosome> new_population;

	for (int i = 0; i < (int)this->population.size(); ++i) {
		
		int r1;
		int r2;
		int r3;

		do {
			r1 = this->randInt((int)this->population.size());
		} while (i == r1);

		do {
			r2 = this->randInt((int)this->population.size());
		} while (i == r2 || r1 == r2);

		do {
			r3 = this->randInt((int)this->population.size());
		} while (i == r3 || r1 == r3 || r2 == r3);

		new_population.push_back(this->mutate(i,r1,r2,r3));
	}

	this->population = new_population;
	this->recordStats();
}

Chromosome EXAFSDE::mutate(int i, int r1, int r2, int r3) {

	std::vector<PDBAtom> modified_atoms;
	Chromosome selectedChromosome = this->population[i];

	int chromosome_length = selectedChromosome.atoms.size();

	// Modify each atom.
	for (int j = 0; j < chromosome_length; ++j) {
		
		PDBAtom temp = this->population[r2].atoms[j] - this->population[r3].atoms[j];
		temp = temp * this->f;
		modified_atoms.push_back(this->population[r1].atoms[j] + temp);
	}

	int Irand = this->randInt(chromosome_length);
	std::vector<PDBAtom> final_atoms;

	// Decide which atoms to replace.
	for (int j = 0; j < chromosome_length; ++j) {
		
		if (this->unifRand() <= this->cr || j == Irand) {
			final_atoms.push_back(modified_atoms[j]);
		} else {
			final_atoms.push_back(selectedChromosome.atoms[j]);
		}
	}

	// Test if fitness has improved.
	Chromosome originalChromosome = selectedChromosome;

	// Replace atoms and evaluate.
	selectedChromosome.atoms = final_atoms;
	this->evaluate(selectedChromosome);

	std::cout << "\t Child: " << i << ", Original: " << originalChromosome.exafs_score << ", Modified: " << selectedChromosome.exafs_score << std::endl;

	if (originalChromosome.exafs_score < selectedChromosome.exafs_score) {
		return originalChromosome;
	} else {
		return selectedChromosome;
	}
}

Chromosome EXAFSDE::best_chromosome() {

	double best = std::numeric_limits<double>::max();
	Chromosome best_chromosome;
	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		if (child->exafs_score < best) {
			best = child->exafs_score;
			best_chromosome = *child;
		}
	}

	return best_chromosome;
}

void EXAFSDE::saveBestChromosome() {

	Chromosome best_individual = this->best_chromosome();
	this->best_individuals.push_back(best_individual);

	std::cout << "Best: " << best_individual.exafs_score << std::endl;
}

double EXAFSDE::unifRand() {
	return rand() / double(RAND_MAX);
}

double EXAFSDE::randInt(int max) {

	return this->unifRand() * max;
}

void EXAFSDE::initStats() {

	this->output_stream.open(( this->stats_folder + "/" + this->results_file).c_str());
	if (this->output_stream.is_open()) {
		std::cout << "Results file ready." << std::endl;
	} else {
		std::cout << "Could not create results file." << std::endl;
	}
}

void EXAFSDE::recordStats() {

	double average_exafs_score = 0;
	double average_potential_energy_score = 0;
	double best_score = std::numeric_limits<double>::max();
	Chromosome best_chromosome;

	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		average_exafs_score += child->exafs_score;
		average_potential_energy_score += child->potential_energy;
		if (child->exafs_score < best_score) {
			best_score = child->exafs_score;
			best_chromosome = *child;
		}
	}

	this->output_stream << best_score << "," << (average_exafs_score/(int)this->population.size()) << "," << best_chromosome.potential_energy << "," << (average_potential_energy_score/(int)this->population.size()) << std::endl;
}

void EXAFSDE::finalStats() {

	this->output_stream.close();

	std::vector< std::pair<double, double> > target_exafs = this->exafs_evaluator->getTargetEXAFS();
	std::ofstream output(this->stats_folder + "/generation_data.csv");
	for (int i = 0; i < (int)this->best_individuals[0].exafs_data.size(); ++i) {
		
		output << this->best_individuals[0].exafs_data[i].first;

		if (this->best_individuals[0].exafs_data[i].first != 0) {
			output << "," << target_exafs[i-1].second;
		} else {
			output << ",0";
		}

		for (int j = 0; j < (int)this->best_individuals.size(); ++j) {
			
			output << "," << this->best_individuals[j].exafs_data[i].second;
		}
		output << std::endl;
	}
	output.close();

	this->exafs_evaluator->updateAtoms(this->best_individuals[this->best_individuals.size()-1].atoms);
	this->exafs_evaluator->writePDB(this->stats_folder + "/best_chromosome.pdb");
	
}

bool EXAFSDE::convergence() {

	std::vector<Chromosome> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_de);
	std::vector<Chromosome>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) == 1;
}

bool EXAFSDE::convergence(double rate) {

	std::vector<Chromosome> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_de);
	std::vector<Chromosome>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) <= (rate * this->population.size());
}