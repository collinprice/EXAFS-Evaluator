
#include "pdbhelper.h"
#include "Genfig/genfig.h"

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

	PDBHelper helper(std::string("relaxed-H.pdb"), std::string("prmtop-sphere"));

	std::cout << helper.numberOfAtoms() << std::endl;

	helper.writePDBFile("tester.txt");

	std::vector<PDBAtom> atoms = helper.getEXAFSAtoms();
	atoms[0].x = 100;
	helper.updateEXAFSAtoms(atoms);
	helper.writePDBFile("tester2.txt");
	
	return 0;
}
