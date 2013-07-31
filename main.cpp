#include "genfig/genfig.h"
#include "iohelper/iohelper.h"
#include "exafsga.h"

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
		srand(time(NULL));
	}

	PDBHelper* pdb_helper = new PDBHelper(config.getString("pdb-file"), config.getString("amber-topology-file"), "temp_pdb.pdb");
	IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(pdb_helper->getEXAFSAtoms(), config.getString("target-atom"), config.getString("experimental-exafs"), config.getFloat("x-min"), config.getFloat("x-max"), config.getString("feff"), config.getString("ifeffit"));
	VMDHelper* vmd_helper = new VMDHelper(pdb_helper->output_pdb_file, config.getString("amber-topology-file"), config.getString("namd2-path"), config.getString("vmd-path"));

	std::vector<std::string> files = IOHelper::contentsOfFolder(config.getString("initial-population-dir"));
	std::vector< std::vector<PDBAtom> > initial_population;

	for (std::vector<std::string>::iterator filename = files.begin(); filename != files.end(); ++filename) {
		
		std::string file = config.getString("initial-population-dir") + "/" + *filename;
		pdb_helper->updateAtomsFromXYZ(file, true);

		initial_population.push_back(pdb_helper->getEXAFSAtoms());
	}

	EXAFSEvaluator* exafs_evaluator = new EXAFSEvaluator(ifeffit_helper, pdb_helper, vmd_helper);

	EXAFSGA ga(exafs_evaluator, config.getFloat("mutation"), config.getFloat("crossover"), config.getBool("elitism"), config.getInt("max-generations"), config.getString("results"));
	ga.begin(initial_population);

	return 0;
}
