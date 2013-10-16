#include "exafsga.h"

#include <iostream>
#include <math.h>
#include <limits>
#include <sstream>

EXAFSGA::EXAFSGA(EXAFSEvaluator* exafs_evaluator, double mutation_rate, double crossover_rate, bool elitism, int max_generations, std::string results_file) {

	this->exafs_evaluator = exafs_evaluator;
	this->mutation_rate = mutation_rate;
	this->crossover_rate = crossover_rate;
	this->elitism = elitism;
	this->max_generations = max_generations;
	this->results_file = results_file;
}

EXAFSGA::~EXAFSGA() {

	delete this->exafs_evaluator;
}

void EXAFSGA::begin(std::vector< std::vector<PDBAtom> > initial_population) {

	this->initPopulation(initial_population);
	this->initStats();
	this->recordStats();

	std::cout << "Begin" << std::endl;
	for (int i = 0; i < this->max_generations; ++i) {
		std::cout << "Generation: " << (i+1) << std::endl;
		this->evolve();

		this->saveBestChromosome();

		if (this->convergence()) break;
	}

	this->finalStats();
}

void EXAFSGA::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Chromosome child(*i);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
	this->saveBestChromosome();
}

void EXAFSGA::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		if (!this->population[i].is_evaluated) {
			this->evaluate(this->population[i]);
			this->population[i].is_evaluated = true;
		}
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << std::endl;
	}
}

void EXAFSGA::evaluate( Chromosome& child ) {

	this->exafs_evaluator->updateAtoms(child.atoms);
	
	child.exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_data = this->exafs_evaluator->getEXAFSData();

	child.potential_energy = this->exafs_evaluator->calculatePotentialEnergy();
}

void EXAFSGA::evolve() {
	
	std::vector<Chromosome> new_population;

	if (this->elitism) {
		// std::cout << "elitism" << std::endl;
		new_population.push_back(this->best_chromosome());
	}

	while ((int)new_population.size() < (int)this->population.size()) {
		
		Chromosome p1 = this->selection();
		Chromosome p2 = this->selection();

		if (this->crossover_rate > (rand() / double(RAND_MAX))) {
			this->crossover(p1,p2);
			p1.is_evaluated = false;
			p2.is_evaluated = false;
		} else if (this->mutation_rate > (rand() / double(RAND_MAX))) {
			this->mutate(p1);
			this->mutate(p2);
			p1.is_evaluated = false;
			p2.is_evaluated = false;
		}

		if ((int)new_population.size() < (int)this->population.size()) {
			new_population.push_back(p1);

			if ((int)new_population.size() != (int)this->population.size()) {
				new_population.push_back(p2);
			}
		}
	}

	this->population = new_population;
	this->evaluatePopulation();
	this->recordStats();
}

Chromosome EXAFSGA::selection() {

	
	// Primitive Tournament Selection
	std::vector<Chromosome> tournament;

	for (int i = 0; i < 3; ++i) {
		int selected = this->population.size() * this->unifRand();
		tournament.push_back(this->population[selected]);
	}

	Chromosome winner;
	double best_score = std::numeric_limits<double>::max();
	for (int i = 0; i < 3; ++i) {
		if (tournament[i].exafs_score < best_score) {
			winner = tournament[i];
			best_score = winner.exafs_score;
		}
	}

	return winner;
}



void EXAFSGA::crossover(Chromosome& p1, Chromosome& p2) {

	// One point crossover
	int pivot = p1.atoms.size() * this->unifRand();

	for (int i = 0; i < pivot; ++i) {
		
		PDBAtom temp = p1.atoms[i];
		p1.atoms[i] = p2.atoms[i];
		p2.atoms[i] = temp;

	}

}

/*	
phi = 0 - 360
theta = 0 - 180

x = r sin theta cos phi
y = r sin theta sin phi
z = r cos theta
*/
void EXAFSGA::mutate(Chromosome& child) {

	int index = child.atoms.size() * this->unifRand();
	double theta = 180 * this->unifRand();
	double phi = 360 * this->unifRand();
	double r = 0.05;

	child.atoms[index].x = r * sin(theta) * cos(phi);
	child.atoms[index].y = r * sin(theta) * sin(phi);
	child.atoms[index].z = r * cos(theta);

}

Chromosome EXAFSGA::best_chromosome() {

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

void EXAFSGA::saveBestChromosome() {

	Chromosome best_individual = this->best_chromosome();
	this->best_individuals.push_back(best_individual);

	std::cout << "Best: " << best_individual.exafs_score << std::endl;
}

double EXAFSGA::unifRand() {
	return rand() / double(RAND_MAX);
}

void EXAFSGA::initStats() {

	this->output_stream.open(this->results_file.c_str());
	if (this->output_stream.is_open()) {
		std::cout << "Results file ready." << std::endl;
	} else {
		std::cout << "Could not create results file." << std::endl;
	}
}

void EXAFSGA::recordStats() {

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

void EXAFSGA::finalStats() {

	this->output_stream.close();

	std::vector< std::pair<double, double> > target_exafs = this->exafs_evaluator->getTargetEXAFS();
	std::ofstream output("generation_data.csv");
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
	this->exafs_evaluator->writePDB("best_chromosome.pdb");
	
}

bool EXAFSGA::convergence() {

	double first_child_score = this->population[0].exafs_score;
	for (int i = 0; i < (int)this->population.size(); ++i) {
		
		if (this->population[i].exafs_score != first_child_score) {

			return false;
		}
	}

	std::cout << "Early convergence. Ending run." << std::endl;
	return true;
}