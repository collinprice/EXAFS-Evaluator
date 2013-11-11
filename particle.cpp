#include "particle.h"

#include <iostream>

// External
double unifRand2(double max) {
	return (rand() / double(RAND_MAX)) * max * 2 - max;
}

double unifRand2() {
	return rand() / double(RAND_MAX);
}

Particle::Particle() {

	this->init();
}

Particle::Particle( std::vector<PDBAtom> atoms, double max_offset ) {

	this->position = atoms;
	this->best_position = atoms;
	this->best_exafs_score = 99999;
	this->init();
	this->generate_initial_velocity(max_offset);
}

Particle::Particle( const Particle& other ) : exafs_score( other.exafs_score ),
													potential_energy( other.potential_energy ),
													is_evaluated( other.is_evaluated ),
													exafs_data( other.exafs_data ),
													position( other.position ),
													velocity( other.velocity ),
													best_exafs_score( other.best_exafs_score ) {}

void Particle::init() {

	this->exafs_score = 0;
	this->potential_energy = 0;
	this->is_evaluated = false;
}

void Particle::generate_initial_velocity(double max_offset) {

	std::vector< std::vector<double> > velocity;
	for (int i = 0; i < (int)this->position.size(); ++i) {
		
		std::vector<double> single;
		single.push_back(unifRand2(max_offset));
		single.push_back(unifRand2(max_offset));
		single.push_back(unifRand2(max_offset));

		velocity.push_back(single);
	}

	this->velocity = velocity;
}

void Particle::update_velocity(std::vector< PDBAtom > global_best_position) {

	double inertia = 2;
	double social = 1.496180;
	double cognitive = 1.496180;

	for (int i = 0; i < (int)this->velocity.size(); ++i) {
		
		this->velocity[i][0] *= (inertia * this->velocity[i][0]) +
			(unifRand2()*social) * global_best_position[i].x - this->position[i].x + 
			(unifRand2()*cognitive) * this->best_position[i].x - this->position[i].x;

		this->velocity[i][1] *= (inertia * this->velocity[i][1]) +
			(unifRand2()*social) * global_best_position[i].y - this->position[i].y + 
			(unifRand2()*cognitive) * this->best_position[i].y - this->position[i].y;

		this->velocity[i][2] *= (inertia * this->velocity[i][2]) +
			(unifRand2()*social) * global_best_position[i].z - this->position[i].z + 
			(unifRand2()*cognitive) * this->best_position[i].z - this->position[i].z;

		std::cout << this->velocity[i][0] << std::endl;
	}
}

void Particle::update_position() {

	for (int i = 0; i < (int)this->position.size(); ++i) {
		
		std::cout << this->velocity[i][0] << std::endl;
		std::cout << this->position[i].x << std::endl;

		this->position[i].x = this->velocity[i][0] + this->position[i].x;
		this->position[i].y = this->velocity[i][1] + this->position[i].y;
		this->position[i].z = this->velocity[i][2] + this->position[i].z;
	}
}

void Particle::update_best_position() {

	if (this->exafs_score < this->best_exafs_score) {
		this->best_exafs_score = this->exafs_score;
		this->best_position = this->position;
	}
}