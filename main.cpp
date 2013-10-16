#include "genfig/genfig.h"
#include "exafsga.h"
#include "dcdhelper.h"
#include "clustering.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>
#include <math.h>

int main(int argc, char **argv) {
	
	int c;
	std::string input_file;
	while ( (c = getopt(argc, argv, "i:")) != -1 ) {

		switch(c) {
			case 'i':
				input_file.assign(optarg);
				break;
			default:
				exit(EXIT_FAILURE);
		}
	}

	if (input_file.size() == 0) {
		std::cout << "Config file required." << std::endl;
		exit(EXIT_FAILURE);
	}

	Genfig config(input_file);


	PDBHelper* pdb_helper = new PDBHelper(config.getString("pdb-file"), config.getString("amber-topology-file"), "temp_pdb.pdb", config.getStringList("exafs-atoms"));
	IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(config.getString("folder-name"), pdb_helper->getEXAFSAtoms(), config.getString("target-atom"), config.getString("experimental-exafs"), config.getDouble("x-min"), config.getDouble("x-max"), config.getString("feff"), config.getString("ifeffit"));
	VMDHelper* vmd_helper = new VMDHelper(pdb_helper->output_pdb_file, config.getString("amber-topology-file"), config.getString("namd2-path"), config.getString("vmd-path"));
	
	double rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);

	pdb_helper->writePDBFile();
	double potential_energy = vmd_helper->calculateEnergy();

	std::cout << "RMSD = " << rmsd << std::endl;
	std::cout << "Potential Energy = " << potential_energy << std::endl;
/*
	EXAFSEvaluator* exafs_evaluator = new EXAFSEvaluator(ifeffit_helper, pdb_helper, vmd_helper);

	EXAFSGA ga(exafs_evaluator, config.getDouble("mutation"), config.getDouble("crossover"), config.getBool("elitism"), config.getInt("max-generations"), config.getString("results"));

	std::cout << "Getting initial population." << std::endl;

	std::vector< std::vector<PDBAtom> > initial_population;
	switch (config.getInt("population-type")) {
		case 0: {
			std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(config.getString("dcd-file"), config.getInt("population-size"));
			for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

				pdb_helper->updateEXAFSAtomsFromXYZ(*i);
				initial_population.push_back( pdb_helper->getEXAFSAtoms() );
			}
			break;
		}
		case 1: {
			std::vector< std::vector<PDBAtom> > initial_exafs_population = random_population(pdb_helper->getEXAFSAtoms(), config.getInt("population-size"), 0.05);
			for (std::vector< std::vector<PDBAtom> >::iterator i = initial_exafs_population.begin(); i != initial_exafs_population.end(); ++i) {
				initial_population.push_back( *i );
			}
			break;
		}
	}

	std::cout << "GA: Begin" << std::endl;
	ga.begin(initial_population);
*/	

	return 0;
}
