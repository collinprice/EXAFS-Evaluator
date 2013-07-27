#include "genfig/genfig.h"
#include "pdbhelper.h"
#include "vmdhelper.h"
#include "ifeffithelper.h"

#include <iostream>
#include <string>
#include <unistd.h>

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

	PDBHelper helper(config.getString("pdb-file"), config.getString("amber-topology-file"));
	std::vector<PDBAtom> original_atoms = helper.getEXAFSAtoms();

	IFEFFITHelper ifeffit_helper(original_atoms, config.getString("target-atom"), config.getString("experimental-exafs"), config.getFloat("x-min"), config.getFloat("x-max"), config.getString("feff"), config.getString("ifeffit"));
	std::cout << ifeffit_helper.run(original_atoms, true) << std::endl;
	ifeffit_helper.clean();

	VMDHelper vmd_helper(config.getString("pdb-file"),config.getString("amber-topology-file"),config.getString("namd2-path"),config.getString("vmd-path"));
	std::cout << vmd_helper.calculateEnergy() << std::endl;

	return 0;
}
