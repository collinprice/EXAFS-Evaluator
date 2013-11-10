#include "genfig/genfig.h"
#include "exafsga.h"
#include "dcdhelper.h"
#include "clustering.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <memory>
#include <math.h>

#include <typeinfo>

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

std::vector<int> dcdGetIndexLessThan(std::string filename, double limit) {

	std::vector<int> indexes;
	std::ifstream index_file(filename.c_str());
	if (index_file.is_open() && index_file.good()) {

		std::string score_str;
		int i = 0;
		while(index_file.good()) {

			index_file >> score_str;
			double score = atof(score_str.c_str());

			if (score < limit) {
				indexes.push_back(i);
			}
			
			++i;
		}
	}

	return indexes;
}

int main(int argc, char **argv) {
	
	int c;
	std::string fitness_file, ga_file;
	while ( (c = getopt(argc, argv, "f:g:")) != -1 ) {

		switch(c) {
			case 'f':
				fitness_file.assign(optarg);
				break;
			case 'g':
				ga_file.assign(optarg);
				break;
			default:
				exit(EXIT_FAILURE);
		}
	}

	if (fitness_file.size() == 0) {
		std::cout << "Fitness file required. Use -f <file>." << std::endl;
		exit(EXIT_FAILURE);
	}

	if (ga_file.size() == 0) {
		std::cout << "GA file required. Use -g <file>." << std::endl;
		exit(EXIT_FAILURE);
	}

	Genfig fitness_config(fitness_file);
	Genfig ga_config(ga_file);

	if (ga_config.hasKey("seed")) {
            srand(ga_config.getInt("seed"));
            std::cout << "Seed: " << ga_config.getInt("seed") << std::endl;
    } else {
            time_t inital_seed = time(NULL);
            srand(inital_seed);
            std::cout << "Seed: " << inital_seed << std::endl;
    }

	PDBHelper* pdb_helper = new PDBHelper(fitness_config.getString("pdb-file"), fitness_config.getString("amber-topology-file"), "temp_pdb.pdb", fitness_config.getStringList("exafs-atoms"));
	IFEFFITHelper* ifeffit_helper = new IFEFFITHelper(fitness_config.getString("folder-name"), pdb_helper->getEXAFSAtoms(), fitness_config.getString("target-atom"), fitness_config.getString("experimental-exafs"), fitness_config.getDouble("x-min"), fitness_config.getDouble("x-max"), fitness_config.getString("feff"), fitness_config.getString("ifeffit"));
	VMDHelper* vmd_helper = new VMDHelper(pdb_helper->output_pdb_file, fitness_config.getString("amber-topology-file"), fitness_config.getString("namd2-path"), fitness_config.getString("vmd-path"));
	EXAFSEvaluator* exafs_evaluator = new EXAFSEvaluator(ifeffit_helper, pdb_helper, vmd_helper);

	if (ga_config.getString("eval-type").compare("solo") == 0) {
		std::cout << "Solo" << std::endl;
		std::cout << "Initial pdb file" << std::endl;
	 	double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
	 	std::cout << "RMSD = " << initial_rmsd << std::endl;

	} else if (ga_config.getString("eval-type").compare("ga") == 0) {
		std::cout << "GA" << std::endl;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::cout << "Getting initial population." << std::endl;

		std::vector< std::vector<PDBAtom> > initial_population;
		switch (ga_config.getInt("population-type")) {
			case 0: {
				std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"), ga_config.getInt("population-size"));
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

					pdb_helper->updateEXAFSAtomsFromXYZ(*i);
					initial_population.push_back( pdb_helper->getEXAFSAtoms() );
				}
				break;
			}
			case 1: {
				std::vector< std::vector<PDBAtom> > initial_exafs_population = random_population(pdb_helper->getEXAFSAtoms(), ga_config.getInt("population-size"), 0.05);
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_exafs_population.begin(); i != initial_exafs_population.end(); ++i) {
					initial_population.push_back( *i );
				}
				break;
			}
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}

		ga.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("ga_recenter") == 0) {
		std::cout << "GA Recentering" << std::endl;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::cout << "Getting initial population." << std::endl;

		std::vector< std::vector<PDBAtom> > initial_population;
		switch (ga_config.getInt("population-type")) {
			case 0: {
				std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"), ga_config.getInt("recentering-population"));
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

					pdb_helper->updateEXAFSAtomsFromXYZ(*i);
					initial_population.push_back( pdb_helper->getEXAFSAtoms() );
				}
				break;
			}
			case 1: {
				std::vector< std::vector<PDBAtom> > initial_exafs_population = random_population(pdb_helper->getEXAFSAtoms(), ga_config.getInt("population-size"), 0.05);
				for (std::vector< std::vector<PDBAtom> >::iterator i = initial_exafs_population.begin(); i != initial_exafs_population.end(); ++i) {
					initial_population.push_back( *i );
				}
				break;
			}
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}

		ga.begin_recentering(initial_populations, ga_config.getInt("population-size"), ga_config.getDouble("convergence-rate"), ga_config.getInt("recentering"));

	} else if (ga_config.getString("eval-type").compare("xyz") == 0) {
		std::cout << "XYZ" << std::endl;

		std::vector<std::string> xyz_files = ga_config.getStringList("xyz-files");

		for (std::vector<std::string>::iterator file = xyz_files.begin(); file != xyz_files.end(); ++file) {
			
			std::cout << "File: " << *file << std::endl;
			std::vector<PDBAtom> atoms = extractAtomsFromXYZ(*file);
			IFEFFITHelper ifeffit_helper = IFEFFITHelper(fitness_config.getString("folder-name"), atoms, fitness_config.getString("target-atom"), fitness_config.getString("experimental-exafs"), fitness_config.getDouble("x-min"), fitness_config.getDouble("x-max"), fitness_config.getString("feff"), fitness_config.getString("ifeffit"));

			double rmsd = ifeffit_helper.run(atoms, true);
			std::cout << "RMSD = " << rmsd << std::endl;

			std::vector< std::pair<double, double> > target_exafs = ifeffit_helper.getTargetEXAFS();
			std::vector< std::pair<double, double> > calculated_exafs = ifeffit_helper.getEXAFSData();

			std::ofstream exafs_output((*file + ".csv").c_str());
			for (int i = 0; i < (int)target_exafs.size() && i < (int)calculated_exafs.size() - 2; ++i) {
				
				exafs_output << target_exafs[i].first << "," << target_exafs[i].second << "," << calculated_exafs[i+1].second << std::endl;
			}
			exafs_output.close();

		}
	} else if (ga_config.getString("eval-type").compare("indexes") == 0) {

		std::cout << "Indexes" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));

	} else if (ga_config.getString("eval-type").compare("index_ga") == 0) {

		std::cout << "Index_GA" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}

		ga.begin(initial_populations);

	} else if (ga_config.getString("eval-type").compare("index_ga_recenter") == 0) {

		std::cout << "Index_GA_Recenter" << std::endl;

		std::vector<int> indexes = dcdGetIndexLessThan(ga_config.getString("index-file"), ga_config.getDouble("index-max"));
		std::vector< std::vector<PDBAtom> > initial_population;

		EXAFSGA ga(exafs_evaluator, ga_config.getDouble("mutation"), ga_config.getDouble("crossover"), ga_config.getBool("elitism"), ga_config.getInt("max-generations"), ga_config.getString("results"));

		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZsByIndex(ga_config.getString("dcd-file"), indexes);
		std::cout << "DCD Population = " << initial_dcd_population.size() << std::endl;
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {

			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			initial_population.push_back( pdb_helper->getEXAFSAtoms() );
		}

		std::cout << "GA: Begin" << std::endl;

		std::vector< std::vector< std::vector<PDBAtom> > > initial_populations;

		for (int i = 0; i < ga_config.getInt("runs"); ++i) {
			initial_populations.push_back(initial_population);
		}
		
		ga.begin_recentering(initial_populations, ga_config.getInt("population-size"), ga_config.getDouble("convergence-rate"), ga_config.getInt("recentering"));

	} else {
		std::cout << "Other" << std::endl;

		// Get all the xyz entries.
		std::vector< std::vector<PDBAtom> > initial_dcd_population = DCDHelper::getXYZs(ga_config.getString("dcd-file"));

		std::ofstream dcd_results((std::string("dcd_results.csv")).c_str());
		for (std::vector< std::vector<PDBAtom> >::iterator i = initial_dcd_population.begin(); i != initial_dcd_population.end(); ++i) {
			pdb_helper->updateEXAFSAtomsFromXYZ(*i);
			double initial_rmsd = ifeffit_helper->run(pdb_helper->getEXAFSAtoms(), true);
			
			std::cout << initial_rmsd << std::endl;
			dcd_results << initial_rmsd << std::endl;
		}
		dcd_results.close();
	}

	return 0;
}
