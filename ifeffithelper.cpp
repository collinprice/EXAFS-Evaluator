#include "ifeffithelper.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <cstdlib>
#include <thread>
#include <unistd.h>
#include <math.h>
#include <limits>
#include <iostream>

const std::string IFEFFITHelper::CALCULATED_EXAFS_FILENAME = "my_chi.chi3";
const std::string IFEFFITHelper::IFEFFIT_SCRIPT = "_ifeffit_script.sh";
const std::string IFEFFITHelper::CLEAN_SCRIPT = "clean.sh";

IFEFFITHelper::IFEFFITHelper(std::vector<PDBAtom> atoms, 
								std::string target_atom, 
								std::string target_exafs_filename,
								double x_min, 
								double x_max,
								std::string feff_path,
								std::string ifeffit_path) {

	std::string feff_header[] = {
		"TITLE My Molecule",
		"POTENTIALS",
		"*	ipot	z	tag"
	};

	std::string feff_mid[] = {
		"CONTROL 1 1 1 1 1 1",
		"NLEG 8",
		"PRINT 0 0 0 3",
		"EXCHANGE 0 0.0 0.0",
		"CRITERIA 4.0 2.5",
		"ATOMS",
		"*	x	y	z	ipot	tag"
	};

	std::string feff_setup_steps[] = {
		feff_path + " > /dev/null",
		"for i in `ls feff00*`",
		"do",
		"\techo \"$i\" >> process_files.dat",
		"done",
		"cd .."
	};

	std::string ifeffit_header[] = {
		"read_data(file=" + target_exafs_filename + ",type=chi,group=data)",
		"my.k   = data.k",
		"my.chi = data.chi/(data.k**3)",
		"newdata.k = range(0,10.0,0.05)",
		"newdata.chi = qinterp(my.k, my.chi, newdata.k)",
		"guess e  =  0.0",
		"set   s  = 1.0",
		"guess s1 = 0.0025"
	};

	std::string ifeffit_footer[] = {
		"set (kmin = 1.0, kmax =10.0)",
		"set (kweight=3,dk = 1, kwindow='hanning')",
		"set (rmin = 0, rmax = 6)",
		"ff2chi(1-100,group = init)",
		"feffit(1-100,chi = newdata.chi, group=fit )",
		"fit3.chi = fit.chi*fit.k**3",
		"write_data(file=my_chi.chi3,fit.k,fit3.chi)",
		"exit"
	};

	std::string ifeffit_script[] = {
		"feff6 > /dev/null",
		ifeffit_path + " -q process.iff"
	};

	int target_atom_atomic_number = PDBAtom::atomicSymbolToNumber(target_atom);

	// Init header array.
	int size = sizeof(feff_header)/sizeof(feff_header[0]);
	for (int i = 0; i < size; ++i) {
		this->cached_header.push_back(feff_header[i]);
	}

	// Find unique atoms and assign indexes.
	int counter = 1;
	for (int i = 0; i < (int)atoms.size(); ++i) {
		if (this->unique_atoms.find(atoms[i].getAtomicNumber()) == unique_atoms.end()) {
			unique_atoms[atoms[i].getAtomicNumber()] = counter++;
		}
		if (atoms[i].getAtomicNumber() == target_atom_atomic_number) {
			target_indexes.push_back(i);
		}
	}

	std::ostringstream oss;
	oss << "\t" << 0 << "\t" << target_atom_atomic_number << "\t" << target_atom;

	this->cached_header.push_back(oss.str());

	for (std::map<int, int>::iterator atom = unique_atoms.begin(); atom != unique_atoms.end(); ++atom) {
		oss.clear();
		oss.str("");

		oss << "\t" << atom->second << "\t" << atom->first << "\t" << PDBAtom::atomicNumberToSymbol(atom->first);
		this->cached_header.push_back(oss.str());
	}

	size = sizeof(feff_mid)/sizeof(feff_mid[0]);
	for (int i = 0; i < size; ++i) {
		this->cached_header.push_back(feff_mid[i]);
	}

	// PREPARE
	std::string setup_feff_script("feff_setup.sh");
	std::ofstream feff_setup(setup_feff_script.c_str());

	// Create folders for each feff/ifeffit run
	
	size = (int)this->target_indexes.size();
	for (int i = 0; i < size; ++i) {
		std::ostringstream oss;
		oss << i;
		mkdir(oss.str().c_str(), 0755);

		feff_setup << "cd " << oss.str().c_str() << std::endl;
		int temp_size = sizeof(feff_setup_steps)/sizeof(feff_setup_steps[0]);
		for (int i = 0; i < temp_size; ++i) {
			feff_setup << feff_setup_steps[i] << std::endl;
		}

		this->cached_feff_filenames.push_back(oss.str() + "/feff.inp");

		std::ofstream ifeffit_script_file((oss.str() + IFEFFITHelper::IFEFFIT_SCRIPT).c_str());
		ifeffit_script_file << "cd " << oss.str() << std::endl;
		temp_size = sizeof(ifeffit_script)/sizeof(ifeffit_script[0]);
		for (int i = 0; i < temp_size; ++i) {
			ifeffit_script_file << ifeffit_script[i] << std::endl;
		}
		ifeffit_script_file.close();

		oss.clear();
		oss.str("");
	}
	feff_setup.close();

	// Create feff files.
	this->updateFEFFFiles(atoms);

	system(("bash " + setup_feff_script).c_str());
	system(("rm " + setup_feff_script).c_str());

	for (int i = 0; i < size; ++i) {
		std::ostringstream oss;
		oss << i;

		std::ifstream process_file((oss.str() + "/process_files.dat").c_str());
		if (!process_file.is_open()) {
			// return false;
		}

		std::vector<std::string> files;
		while(process_file.good()) {
			std::string name;
			process_file >> name;
			if (name.size() == 0) continue;
			files.push_back(name);
		}
		process_file.close();
		unlink((oss.str() + "/process_files.dat").c_str());


		std::ofstream ifeffit_file((oss.str() + "/process.iff").c_str());
		if (!ifeffit_file.is_open()) {
			// return false;
		}

		int temp_size = sizeof(ifeffit_header)/sizeof(ifeffit_header[0]);
		for (int i = 0; i < temp_size; ++i) {
			ifeffit_file << ifeffit_header[i] << std::endl;
		}

		for (int i = 0; i < (int)files.size(); ++i) {
			ifeffit_file << "path(" << (i+1) << ", file=" << files[i] << ", e0 = e, s02 = s, sigma2 = s1)" << std::endl;
		}

		temp_size = sizeof(ifeffit_footer)/sizeof(ifeffit_footer[0]);
		for (int i = 0; i < temp_size; ++i) {
			ifeffit_file << ifeffit_footer[i] << std::endl;
		}

		ifeffit_file.close();

		// Move target exafs file into folder.
		system(("cp " + target_exafs_filename + " " + oss.str()).c_str());

		oss.clear();
		oss.str("");
	}

	this->readTargetEXAFS(target_exafs_filename);
	this->x_min = x_min;
	this->x_max = x_max;

	this->clean_script();
}

IFEFFITHelper::~IFEFFITHelper() {
	this->clean();
}

double IFEFFITHelper::run(std::vector<PDBAtom> updated_atoms, bool threaded) {

	this->removeAllCalculatedEXAFSFiles(); // Clean up before run.
	this->updateFEFFFiles(updated_atoms);
	this->processIFEFFIT(threaded);

	return this->calculateRMSD();
}

std::vector< std::pair<double, double> > IFEFFITHelper::getEXAFSData() {

	return this->averaged_calculated_data;
}

std::vector< std::pair<double, double> > IFEFFITHelper::getTargetEXAFS() {
	return this->cached_target_exafs;
}

void IFEFFITHelper::clean() {
	system(("bash " + IFEFFITHelper::CLEAN_SCRIPT).c_str());
}

// PRIVATE FUNCTIONS

void IFEFFITHelper::readTargetEXAFS(std::string filename) {

	std::ifstream exafs_file(filename.c_str());

	if (exafs_file.is_open() && exafs_file.good()) {

		std::string x,y;
		while(exafs_file.good()) {
			exafs_file >> x >> y;
			cached_target_exafs.push_back(std::make_pair(atof(x.c_str()), atof(y.c_str())));
		}

		exafs_file.close();
	}
}

bool IFEFFITHelper::generateFEFFFile(std::vector<PDBAtom> atomic_coordinates, int index, std::string filename) {

	std::string feff_footer[] = {
		"END"
	};

	std::ofstream feff_file(filename.c_str());
	if (feff_file.is_open() && feff_file.good()) {
		
		// Write out header
		for (std::vector<std::string>::iterator i = this->cached_header.begin(); i != this->cached_header.end(); ++i) {
			feff_file << *i << std::endl;
		}	

		int size = (int)atomic_coordinates.size();
		for (int i = 0; i < size; ++i) {

			feff_file << std::setprecision(15) << "\t" << atomic_coordinates[i].x << "\t" << atomic_coordinates[i].y << "\t" << atomic_coordinates[i].z << "\t";

			if (i == index) {
				feff_file << 0;
			} else {
				feff_file << this->unique_atoms[atomic_coordinates[i].getAtomicNumber()];
			}

			feff_file << "\t" << atomic_coordinates[i].atomic_symbol << std::endl;
		}

		size = sizeof(feff_footer)/sizeof(feff_footer[0]);
		for (int i = 0; i < size; ++i) {
			feff_file << feff_footer[i] << std::endl;
		}

		feff_file.close();
		return true;
	} else {
		return false;
	}
}

bool IFEFFITHelper::updateFEFFFiles(std::vector<PDBAtom> atomic_coordinates) {

	int size = (int)this->target_indexes.size();
	for (int i = 0; i < size; ++i) {
		if(!this->generateFEFFFile(atomic_coordinates, this->target_indexes[i], this->cached_feff_filenames[i])) return false;
	}

	return true;
}

void IFEFFITHelper::processIFEFFIT(bool threaded) {

	if (threaded) {
		std::vector<std::string> commands;
		std::vector<std::thread> threads;
		for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
			std::ostringstream oss;
			oss << i;
			std::string bash_string = "bash " + oss.str() + IFEFFITHelper::IFEFFIT_SCRIPT;
			commands.push_back(bash_string);
		}

		for (int i = 0; i < (int)commands.size(); ++i) {
			threads.push_back(std::thread(IFEFFITHelper::staticEntry, commands[i].c_str()));
		}

		for (int i = 0; i < (int)threads.size(); ++i) {
			threads[i].join();
		}
	} else {
		
		for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
			std::ostringstream oss;
			oss << i;
			system(("bash " + oss.str() + IFEFFITHelper::IFEFFIT_SCRIPT).c_str());
		}
	}
	
}

void IFEFFITHelper::staticEntry(const char* command) {
	system(command);
}

double IFEFFITHelper::calculateRMSD() {

	// Check that all output files are there. If there are any missing then ifeffit failed.
	if (!this->canPerformIFEFFITCalculations()) {
		return std::numeric_limits<double>::max();
	}

	// Read in target_indexes output files.
	// Average all files.
	std::vector< std::vector< std::pair<double, double> > > calculated_exafs;
	std::ostringstream oss;

	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		std::vector< std::pair<double, double> > exafs_data;
		oss << i;
		std::ifstream exafs_file(oss.str() + "/" + IFEFFITHelper::CALCULATED_EXAFS_FILENAME);

		std::string x,y;
		// Skip the first three lines. They are headers.
		std::getline(exafs_file, x);
		std::getline(exafs_file, x);
		std::getline(exafs_file, x);
		while(exafs_file.good()) {
			exafs_file >> x >> y;
			exafs_data.push_back( std::make_pair(atof(x.c_str()), atof(y.c_str())) );
		}
		calculated_exafs.push_back(exafs_data);

		exafs_file.close();
		oss.clear();
		oss.str("");
	}

	this->averaged_calculated_data.clear();

	// Preload with 0's
	for (int i = 0; i < (int)calculated_exafs[0].size(); ++i) {
		this->averaged_calculated_data.push_back( std::make_pair(calculated_exafs[0][i].first, 0) );
	}

	for (int i = 0; i < (int)calculated_exafs.size(); ++i) {
		for (int j = 0; j < (int)calculated_exafs[i].size(); ++j) {
			this->averaged_calculated_data[j].second += calculated_exafs[i][j].second;
		}
	}

	for (int i = 0; i < (int)this->averaged_calculated_data.size(); ++i) {
		this->averaged_calculated_data[i].second /= calculated_exafs.size();
	}

	std::vector< std::pair<double, double> > calc_data = this->averaged_calculated_data;
	std::vector< std::pair<double, double> > exp_data = this->cached_target_exafs;

	double rmsd = 0;
	int calc_index = 0;
	int exp_index = 0;
	while (calc_index < (int)calc_data.size() && exp_index < (int)exp_data.size()) {

		if (fabs(calc_data[calc_index].first - exp_data[exp_index].first) < std::numeric_limits<double>::epsilon()) {
			if (calc_data[calc_index].first > this->x_min && calc_data[calc_index].first < this->x_max) {
				double scale_calc = calc_data[calc_index].second;
				double scale_exp = exp_data[exp_index].second;
				rmsd += pow(scale_calc - scale_exp,2);
			}
			++calc_index;
			++exp_index;
		} else if (calc_data[calc_index].first < exp_data[exp_index].first) {
			++calc_index;
		} else {
			++exp_index;
		}
	}

	return rmsd;
}

void IFEFFITHelper::clean_script() {

	std::ofstream cleaner(IFEFFITHelper::CLEAN_SCRIPT);
	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		cleaner << "rm -rf " << i << std::endl;
		cleaner << "rm -rf " << i << IFEFFITHelper::IFEFFIT_SCRIPT << std::endl;
	}
	cleaner << "rm " << IFEFFITHelper::CLEAN_SCRIPT << std::endl;
	cleaner.close();
}

bool IFEFFITHelper::canPerformIFEFFITCalculations() {

	bool has_files = true;
	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		std::ostringstream oss;
		oss << i;
		std::string filename = oss.str() + "/" + IFEFFITHelper::CALCULATED_EXAFS_FILENAME;
		
		std::ifstream temp_file(filename.c_str());
		if (!temp_file.good()) {
			has_files = false;
			std::cout << "ERROR" << std::endl;
		} else {
			temp_file.close();
		}
	}

	return has_files;
}

void IFEFFITHelper::removeAllCalculatedEXAFSFiles() {

	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		std::ostringstream oss;
		oss << i;
		std::string filename = oss.str() + "/" + IFEFFITHelper::CALCULATED_EXAFS_FILENAME;
		remove(filename.c_str());
	}
}