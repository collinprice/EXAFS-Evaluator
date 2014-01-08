#pragma once
#include "ifeffithelper.h"
#include "pdbhelper.h"
#include "vmdhelper.h"

class EXAFSEvaluator {

public:

	EXAFSEvaluator(IFEFFITHelper* ifeffit_helper, PDBHelper* pdb_helper, VMDHelper* vmd_helper);
	~EXAFSEvaluator();

	std::vector<PDBAtom> getAtoms();
	void updateAtoms(std::vector<PDBAtom> atoms);
	double calculateRMSD();
	double calculatePotentialEnergy();
	std::vector< std::pair<double, double> > getEXAFSData();
	std::vector< std::pair<double, double> > getTargetEXAFS();
	void writePDB(std::string filename);

private:

	IFEFFITHelper* ifeffit_helper;
	PDBHelper* pdb_helper;
	VMDHelper* vmd_helper;

};