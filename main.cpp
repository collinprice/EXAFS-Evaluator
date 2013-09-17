#include "genfig/genfig.h"
#include "exafsga.h"
#include "dcdhelper.h"
#include "clustering.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>

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

	if (config.hasKey("seed")) {
		srand(config.getInt("seed"));
	} else {
		time_t inital_seed = time(NULL);
		srand(inital_seed);
		std::cout << "Seed: " << inital_seed << std::endl;
	}

	PDBHelper* pdb_helper = new PDBHelper(config.getString("pdb-file"), config.getString("amber-topology-file"), "temp_pdb.pdb", config.getStringList("exafs-atoms"));
	IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(config.getString("folder-name"), pdb_helper->getEXAFSAtoms(), config.getString("target-atom"), config.getString("experimental-exafs"), config.getDouble("x-min"), config.getDouble("x-max"), config.getString("feff"), config.getString("ifeffit"));
	VMDHelper* vmd_helper = new VMDHelper(pdb_helper->output_pdb_file, config.getString("amber-topology-file"), config.getString("namd2-path"), config.getString("vmd-path"));
	EXAFSEvaluator* exafs_evaluator = new EXAFSEvaluator(ifeffit_helper, pdb_helper, vmd_helper);

	EXAFSGA ga(exafs_evaluator, config.getDouble("mutation"), config.getDouble("crossover"), config.getBool("elitism"), config.getInt("max-generations"), config.getString("results"));

	std::cout << "Getting inital population." << std::endl;
	std::vector< std::vector<PDBAtom> > initial_population = DCDHelper::getXYZs(config.getString("dcd-file"), config.getDouble("percent-of-dcd"));

	// std::vector< std::vector<PDBAtom> > updatedEXAFSPopulation;
	// for (std::vector< std::vector<PDBAtom> >::iterator ind = initial_population.begin(); ind != initial_population.end(); ++ind) {
		
	// 	pdb_helper->updateAtomsFromList(*ind);
	// 	updatedEXAFSPopulation.push_back(pdb_helper->getEXAFSAtoms());
	// }

	// std::vector< std::vector<double> > table = Clustering::createTable(updatedEXAFSPopulation);

	// std::ofstream output("table.csv");
	// for (int i = 0; i < (int)table.size(); ++i) {
	// 	for (int j = 0; j < (int)table.size(); ++j) {
			
	// 		output << table[i][j] << ",";
	// 		if (i == j) break;
	// 	}
	// 	output << std::endl;
	// }
	// output.close();

	std::cout << "GA: Begin" << std::endl;
	ga.begin(initial_population);

	return 0;
}
