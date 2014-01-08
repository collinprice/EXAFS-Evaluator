#include "chromosome.h"

Chromosome::Chromosome() {

	this->init();
}

Chromosome::Chromosome( std::vector<PDBAtom> atoms ) {

	this->atoms = atoms;
	this->init();
}

Chromosome::Chromosome( const Chromosome& other ) : exafs_score( other.exafs_score ),
													potential_energy( other.potential_energy ),
													is_evaluated( other.is_evaluated ),
													exafs_data( other.exafs_data ),
													atoms( other.atoms ) {}

void Chromosome::init() {

	this->exafs_score = 0;
	this->potential_energy = 0;
	this->is_evaluated = false;
}
