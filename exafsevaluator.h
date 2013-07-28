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

private:

	IFEFFITHelper* ifeffit_helper;
	PDBHelper* pdb_helper;
	VMDHelper* vmd_helper;

};