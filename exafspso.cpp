#include "exafspso.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>
#include <sys/stat.h>
#include <algorithm>

EXAFSPSO::EXAFSPSO(EXAFSEvaluator* exafs_evaluator, double inertia, double social, double cognitive, double velocity_range, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->velocity_range = velocity_range;
	this->inertia = inertia;
	this->social = social;
	this->cognitive = cognitive;
	this->max_generations = max_generations;
	this->results_file = results_file;
}

EXAFSPSO::~EXAFSPSO() {

	delete this->exafs_evaluator;
}

bool chromosome_sort_pso(Particle const & a, Particle const & b) {
	return a.exafs_score < b.exafs_score;
}

void EXAFSPSO::begin(std::vector< std::vector< std::vector<PDBAtom> > > initial_populations) {

	for (int i = 0; i < (int)initial_populations.size(); ++i) {
		
		std::stringstream ss;
		ss << (i+1);
		this->stats_folder = "run" + ss.str();
		mkdir(this->stats_folder.c_str(), 0755);

		this->initPopulation(initial_populations[i]);
		this->initStats();
		this->recordStats();

		this->global_best_particle = this->best_particle();
		this->updateGlobalBest();

		std::cout << "Begin Run " << (i+1) << std::endl;
		for (int i = 0; i < this->max_generations; ++i) {
			std::cout << "Generation: " << (i+1) << std::endl;

			this->updateVelocities();
			this->updatePositions();

			this->evaluatePopulation();
			this->recordStats();

			this->updateGlobalBest();

			this->saveBestParticle();
			if (this->convergence()) break;
		}

		this->finalStats();
	}
}

void EXAFSPSO::updateVelocities() {

	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->population[i].updateVelocity(this->global_best_particle,this->inertia,this->social,this->cognitive);
	}
}

void EXAFSPSO::updatePositions() {

	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->population[i].updatePosition();
	}
}

void EXAFSPSO::updateGlobalBest() {

	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->population[i].updateBest();
		if (this->population[i].exafs_score < this->global_best_particle.exafs_score) {
			this->global_best_particle = this->population[i];
		}
	}
}

void EXAFSPSO::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	this->best_individuals.clear();
	this->population.clear();
	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Particle child(this->velocity_range,*i);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
	this->saveBestParticle();
}

void EXAFSPSO::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		if (!this->population[i].is_evaluated) {
			this->evaluate(this->population[i]);
			this->population[i].is_evaluated = true;
		}
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << std::endl;
	}
}

void EXAFSPSO::evaluate( Particle& child ) {

	this->exafs_evaluator->updateAtoms(child.atoms);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

Particle EXAFSPSO::best_particle() {

	double best = std::numeric_limits<double>::max();
	Particle best_particle;
	for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		if (child->exafs_score < best) {
			best = child->exafs_score;
			best_particle = *child;
		}
	}

	return best_particle;
}

void EXAFSPSO::saveBestParticle() {

	Particle best_individual = this->best_particle();
	this->best_individuals.push_back(best_individual);

	std::cout << "Best: " << best_individual.exafs_score << std::endl;
	std::cout << "Global Best: " << global_best_particle.exafs_score << std::endl;
}

double EXAFSPSO::unifRand() {
	return rand() / double(RAND_MAX);
}

double EXAFSPSO::randInt(int max) {

	return this->unifRand() * max;
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
	double best_score = std::numeric_limits<double>::max();
	Particle best_chromosome;

	for (std::vector<Particle>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		average_exafs_score += child->exafs_score;
		if (child->exafs_score < best_score) {
			best_score = child->exafs_score;
			best_chromosome = *child;
		}
	}

	this->output_stream << best_score << "," << (average_exafs_score/(int)this->population.size()) << "," << best_chromosome.potential_energy << std::endl;
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
			
			if ((int)this->best_individuals[j].exafs_data.size() == j) continue; // guard for no exafs data. occurs when no individuals are valid.

			output << "," << this->best_individuals[j].exafs_data[i].second;
		}
		output << std::endl;
	}
	output.close();

	this->exafs_evaluator->updateAtoms(this->best_individuals[this->best_individuals.size()-1].atoms);
	this->exafs_evaluator->writePDB(this->stats_folder + "/best_chromosome.pdb");
	
}

bool EXAFSPSO::convergence() {

	std::vector<Particle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_pso);
	std::vector<Particle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) == 1;
}

bool EXAFSPSO::convergence(double rate) {

	std::vector<Particle> population_copy = this->population;
	std::sort(population_copy.begin(), population_copy.end(), chromosome_sort_pso);
	std::vector<Particle>::iterator it = std::unique(population_copy.begin(), population_copy.end());

	return std::distance(population_copy.begin(), it) <= (rate * this->population.size());
}
