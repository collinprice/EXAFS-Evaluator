#include "particle.h"
#include <iostream>

Particle::Particle() : Chromosome() {

	this->best_exafs_score = 0;
	this->best_potential_energy = 0;
}

Particle::Particle(double range, std::vector<PDBAtom> atoms) : Chromosome(atoms) {

	for (int i = 0; i < (int)atoms.size(); ++i) {
		velocity.push_back(PDBAtom(this->randomRange(range),this->randomRange(range),this->randomRange(range)));
	}

	this->best_exafs_score = -1;
	this->best_potential_energy = -1;
}

Particle::Particle(const Particle& other) 
: Chromosome(other),
  best_exafs_score( other.best_exafs_score ),
  best_potential_energy( other.best_potential_energy ),
  best_exafs_data( other.best_exafs_data ),
  best_atoms( other.best_atoms ),
  velocity( other.velocity ) {}

void Particle::updateBest() {

	std::cout << "updating personal best" << std::endl;
	if (this->exafs_score < this->best_exafs_score) {
		std::cout << "updated" << std::endl;
		this->best_exafs_score = this->exafs_score;
		this->best_exafs_data = this->exafs_data;
		this->best_atoms = this->atoms;
	}
}

void Particle::updatePosition() {

	for (int i = 0; i < (int)this->atoms.size(); ++i) {
		this->atoms[i] = this->atoms[i] + this->velocity[i];
	}
}

void Particle::updateVelocity(Particle global_best) {

	double inertia = 2;
	double social = 1.496180;
	double cognitive = social;
	
	for (int i = 0; i < (int)this->atoms.size(); ++i) {
		
		std::cout << "1: " << this->best_atoms[i].x << std::endl;
		std::cout << "2: " << this->atoms[i].x << std::endl;

		this->velocity[i] = this->velocity[i]*inertia 
		+ (this->best_atoms[i]-this->atoms[i])*social*this->unifRand()
		+ (global_best.atoms[i]-this->atoms[i])*cognitive*this->unifRand();
	}
}

double Particle::unifRand() {
	return rand() / double(RAND_MAX);
}

double Particle::randomRange(double max) {

	double temp = max*this->unifRand();
	return this->unifRand() > 0.5 ? temp : -temp;
}