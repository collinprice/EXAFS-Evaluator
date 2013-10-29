#include "genfig/genfig.h"
#include "exafsga.h"
#include "dcdhelper.h"
#include "clustering.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>
#include <math.h>

std::vector<PDBAtom> extractAtomsFromXYZ(std::string filename) {

	std::vector<PDBAtom> atoms;
	std::ifstream xyz_file(filename.c_str());
	if (xyz_file.is_open() && xyz_file.good()) {

		std::string atomic_symbol, x, y, z;
		while(xyz_file.good()) {

			xyz_file >> atomic_symbol >> x >> y >> z;
			atoms.push_back(PDBAtom(atomic_symbol, atof(x.c_str()), atof(y.c_str()), atof(z.c_str())));
		}

	}
	xyz_file.close();

	return atoms;
}

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
	// VMDHelper* vmd_helper = new VMDHelper(pdb_helper->output_pdb_file, config.getString("amber-topology-file"), config.getString("namd2-path"), config.getString("vmd-path"));
	
	std::cout << "Initial pdb file" << std::endl;
 	double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
 	std::cout << "RMSD = " << initial_rmsd << std::endl;

	// std::vector<std::string> xyz_files = config.getStringList("xyz-files");

	// for (std::vector<std::string>::iterator file = xyz_files.begin(); file != xyz_files.end(); ++file) {
		
	// 	std::cout << "File: " << *file << std::endl;
	// 	std::vector<PDBAtom> atoms = extractAtomsFromXYZ(*file);
	// 	IFEFFITHelper ifeffit_helper = IFEFFITHelper(config.getString("folder-name"), atoms, config.getString("target-atom"), config.getString("experimental-exafs"), config.getDouble("x-min"), config.getDouble("x-max"), config.getString("feff"), config.getString("ifeffit"));

	// 	double rmsd = ifeffit_helper.run(atoms, true);
	// 	std::cout << "RMSD = " << rmsd << std::endl;

	// }

	return 0;
}
