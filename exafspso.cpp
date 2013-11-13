#include "exafspso.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

EXAFSPSO::EXAFSPSO(EXAFSEvaluator* exafs_evaluator, double mutation_rate, double crossover_rate, bool elitism, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->mutation_rate = mutation_rate;
	this->crossover_rate = crossover_rate;
	this->elitism = elitism;
	this->max_generations = max_generations;
	this->results_file = results_file;
}

EXAFSPSO::~EXAFSPSO() {

	delete this->exafs_evaluator;
}

void EXAFSPSO::begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations, double max_initial_velocty_offset) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {
		
		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		this->initPopulation(initial_populations[i], max_initial_velocty_offset);
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

void EXAFSPSO::initPopulation(std::vector< std::vector<PDBAtom> > population, double max_initial_velocty_offset) {

	this->best_individuals.clear();
	this->population.clear();
	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Particle child(*i, max_initial_velocty_offset);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
	this->saveBestChromosome();
}

void EXAFSPSO::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->evaluate(this->population[i]);
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << std::endl;
	}
}

void EXAFSPSO::evaluate( Particle& child ) {

	this->exafs_evaluator->updateAtoms(child.position);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	child.update_best_position();
	// child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

void EXAFSPSO::evolve() {

	for (std::vector<Particle>::iterator individual = this->population.begin(); individual != this->population.end(); ++individual) {
		
		individual->update_velocity(this->best_particle.position);
		individual->update_position();
	}

	this->evaluatePopulation();
	this->recordStats();
}

Particle EXAFSPSO::best_chromosome() {

	double best = std::numeric_limits<double>::max();
	Particle best_chromosome;
	for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		if (child->exafs_score < best) {
			best = child->exafs_score;
			best_chromosome = *child;
		}
	}

	return best_chromosome;
}

void EXAFSPSO::saveBestChromosome() {

	Particle best_individual = this->best_chromosome();
	this->best_individuals.push_back(best_individual);

	this->best_particle = best_individual;

	std::cout << "Best: " << best_individual.exafs_score << std::endl;
}

double EXAFSPSO::unifRand() {
	return rand() / double(RAND_MAX);
}

void EXAFSPSO::initStats() {

	this->output_stream.open(( this->stats_folder + "/" + this->results_file).c_str());
	if (this->output_stream.is_open()) {
		std::cout << "Results file ready." << std::endl;
	} else {
		std::cout << "Could not create results file." << std::endl;
	}
}

void EXAFSPSO::recordStats() {

	double average_exafs_score = 0;
	double average_potential_energy_score = 0;
	double best_score = std::numeric_limits<double>::max();
	Particle best_chromosome;

	for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		average_exafs_score += child->exafs_score;
		average_potential_energy_score += child->potential_energy;
		if (child->exafs_score < best_score) {
			best_score = child->exafs_score;
			best_chromosome = *child;
		}
	}

	this->output_stream << best_score << "," << (average_exafs_score/(int)this->population.size()) << "," << best_chromosome.potential_energy << "," << (average_potential_energy_score/(int)this->population.size()) << std::endl;
}

void EXAFSPSO::finalStats() {

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

	this->exafs_evaluator->updateAtoms(this->best_individuals[this->best_individuals.size()-1].position);
	this->exafs_evaluator->writePDB(this->stats_folder + "/best_chromosome.pdb");
	
}

bool chromosome_sort(Particle const & a, Particle const & b) {
	return a.exafs_score < b.exafs_score;
}

bool EXAFSPSO::convergence() {

	std::vector<Particle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort);
	std::vector<Particle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) == 1;
}

bool EXAFSPSO::convergence(double rate) {

	std::vector<Particle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort);
	std::vector<Particle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) <= (rate * this->population.size());
}