#include "exafsga.h"

#include <iostream>
#include <math.h>
#include <limits>

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
		std::cout << "Best: " << this->best_chromosome().exafs_score << std::endl;
	}

	this->closeStats();
}

void EXAFSGA::initPopulation(std::vector< std::vector<PDBAtom> > population) {

	for (std::vector< std::vector<PDBAtom> >::iterator i = population.begin(); i != population.end(); ++i) {
		
		Chromosome child(*i);
		this->population.push_back(child);
	}

	this->evaluatePopulation();
}

void EXAFSGA::evaluatePopulation() {
	for (int i = 0; i < (int)this->population.size(); ++i) {
		this->evaluate(this->population[i]);
		std::cout << "\t Child: " << i << ", " << this->population[i].exafs_score << std::endl;
	}
}

void EXAFSGA::evaluate( Chromosome& child ) {
	this->exafs_evaluator->updateAtoms(child);
	
	double exafs_score = this->exafs_evaluator->calculateRMSD();
	child.exafs_score = exafs_score;
}

void EXAFSGA::evolve() {
	
	std::vector<Chromosome> new_population;

	if (this->elitism) {
		new_population.push_back(this->best_chromosome());
	}

	while ((int)new_population.size() < (int)this->population.size()) {
		
		Chromosome p1 = this->selection();
		Chromosome p2 = this->selection();

		if (this->crossover_rate > (rand() / double(RAND_MAX))) {
			this->crossover(p1,p2);
		} else if (this->mutation_rate > (rand() / double(RAND_MAX))) {
			this->mutate(p1);
			this->mutate(p2);
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
	int pivot = p1.size() * this->unifRand();

	for (int i = 0; i < pivot; ++i) {
		
		PDBAtom temp = p1[i];
		p1[i] = p2[i];
		p2[i] = temp;

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

	int index = child.size() * this->unifRand();
	double theta = 180 * this->unifRand();
	double phi = 360 * this->unifRand();
	double r = 0.05;

	child[index].x = r * sin(theta) * cos(phi);
	child[index].y = r * sin(theta) * sin(phi);
	child[index].z = r * cos(theta);

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

	double average_score = 0;
	double best_score = std::numeric_limits<double>::max();

	for (std::vector<Chromosome>::iterator child = this->population.begin(); child != this->population.end(); ++child) {
		
		average_score += child->exafs_score;
		if (child->exafs_score < best_score) {
			best_score = child->exafs_score;
		}
	}

	this->output_stream << best_score << "," << (average_score/(int)this->population.size()) << std::endl;
}

void EXAFSGA::closeStats() {

	this->output_stream.close();
}