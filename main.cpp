#include "genfig/genfig.h"
#include "exafsga.h"
#include "dcdhelper.h"
#include "clustering.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>
#include <math.h>

double unifRand() {
        return rand() / double(RAND_MAX);
}

PDBAtom random_move(PDBAtom atom, double max_radius) {

        double theta = 180 * unifRand();
        double phi = 360 * unifRand();
        double r = max_radius * unifRand();

        atom.x = r * sin(theta) * cos(phi);
        atom.y = r * sin(theta) * sin(phi);
        atom.z = r * cos(theta);

        return atom;
}

std::vector< std::vector<PDBAtom> > random_population(std::vector<PDBAtom> seed, int size, double max_radius) {

        std::cout << "Size = " << size << std::endl;

        std::vector< std::vector<PDBAtom> > population;

        for (int i = 0; i < size; ++i) {
                
                std::vector<PDBAtom> individual = seed;

                for (int i = 0; i < (int)individual.size(); ++i) {
                        individual[i] = random_move(individual[i], max_radius);
                }

                population.push_back(individual);
        }

        return population;
}

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

	if (config.hasKey("seed")) {
            srand(config.getInt("seed"));
            std::cout << "Seed: " << config.getInt("seed") << std::endl;
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

	// std::cout << "Initial pdb file" << std::endl;
 // 	double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
 // 	std::cout << "RMSD = " << initial_rmsd << std::endl;

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
