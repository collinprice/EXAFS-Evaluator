#include "exafsevaluator.h"

EXAFSEvaluator::EXAFSEvaluator(IFEFFITHelper* ifeffit_helper, PDBHelper* pdb_helper, VMDHelper* vmd_helper) {

	this->ifeffit_helper = ifeffit_helper;
	this->pdb_helper = pdb_helper;
	this->vmd_helper = vmd_helper;
}

EXAFSEvaluator::~EXAFSEvaluator() {

	delete this->ifeffit_helper;
	delete this->pdb_helper;
	delete this->vmd_helper;
}

std::vector<PDBAtom> EXAFSEvaluator::getAtoms() {

	return this->pdb_helper->getEXAFSAtoms();
}

void EXAFSEvaluator::updateAtoms(std::vector<PDBAtom> atoms) {

	this->pdb_helper->updateEXAFSAtoms(atoms);
}

double EXAFSEvaluator::calculateRMSD() {

	return this->ifeffit_helper->run(this->getAtoms(), true);
}

double EXAFSEvaluator::calculatePotentialEnergy() {

	return this->vmd_helper->calculateEnergy();
}
