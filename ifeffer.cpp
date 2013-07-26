#include "ifeffer.h"
#include "point.h"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <sys/stat.h>
#include <cstdlib>
#include <thread>
#include <unistd.h>
#include <math.h>

#ifdef DEBUG
#include <iostream>
#endif

const std::string Ifeffer::IFEFFIT_SCRIPT = "_ifeffit_script.sh";

Ifeffer::Ifeffer() {

}

Ifeffer::Ifeffer(std::vector<Atom> atoms, Atom target_atom) {

	this->atoms = atoms;

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

	int size = sizeof(feff_header)/sizeof(feff_header[0]);
	for (int i = 0; i < size; ++i) {
		this->cached_header.push_back(feff_header[i]);
	}

	// Find unique atoms and assign indexes.
	int counter = 1;
	for (int i = 0; i < (int)atoms.size(); ++i) {
		if (this->unique_atoms.find(atoms[i].atomic_number) == unique_atoms.end()) {
			unique_atoms[atoms[i].atomic_number] = counter++;
		}
		if (atoms[i].atomic_number == target_atom.atomic_number) {
			target_indexes.push_back(i);
		}
	}

	std::ostringstream oss;
	oss << "\t" << 0 << "\t" << target_atom.atomic_number << "\t" << target_atom.atomic_symbol;

	this->cached_header.push_back(oss.str());

	for (std::map<int, int>::iterator atom = unique_atoms.begin(); atom != unique_atoms.end(); ++atom) {
		oss.clear();
		oss.str("");

		Atom temp_atom = Atom(atom->first);

		oss << "\t" << atom->second << "\t" << temp_atom.atomic_number << "\t" << temp_atom.atomic_symbol;
		this->cached_header.push_back(oss.str());
	}

	size = sizeof(feff_mid)/sizeof(feff_mid[0]);
	for (int i = 0; i < size; ++i) {
		this->cached_header.push_back(feff_mid[i]);
	}

	#ifdef DEBUG
	std::cout << "FEFF Header" << std::endl;
	for (std::vector<std::string>::iterator i = this->cached_header.begin(); i != this->cached_header.end(); ++i) {
		std::cout << *i << std::endl;
	}
	#endif
}

void Ifeffer::readTargetEXAFS(std::string filename) {

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

bool Ifeffer::generateFEFFFile(std::vector<Point> atomic_coordinates, int index, std::string filename) {

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
			Atom temp_atom = this->atoms[i];

			feff_file << std::setprecision(15) << "\t" << atomic_coordinates[i].x << "\t" << atomic_coordinates[i].y << "\t" << atomic_coordinates[i].z << "\t";

			if (i == index) {
				feff_file << 0;
			} else {
				feff_file << this->unique_atoms[temp_atom.atomic_number];
			}

			feff_file << "\t" << temp_atom.atomic_symbol << std::endl;
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

bool Ifeffer::updateFEFFFiles(std::vector<Point> atomic_coordinates) {

	int size = (int)this->target_indexes.size();
	for (int i = 0; i < size; ++i) {
		if(!this->generateFEFFFile(atomic_coordinates, this->target_indexes[i], this->cached_feff_filenames[i])) return false;
	}

	return true;
}

void Ifeffer::processIFEFFIT(bool threaded) {

	if (threaded) {
		std::vector<std::string> commands;
		std::vector<std::thread> threads;
		for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
			std::ostringstream oss;
			oss << i;

			std::string bash_string = "bash " + oss.str() + Ifeffer::IFEFFIT_SCRIPT;
			commands.push_back(bash_string);
			#ifdef DEBUG
			std::cout << "Creating: " << bash_string << std::endl;
			#endif

			oss.clear();
			oss.str("");
		}

		for (int i = 0; i < (int)commands.size(); ++i) {
			threads.push_back(std::thread(Ifeffer::staticEntry, commands[i].c_str()));
		}

		for (int i = 0; i < (int)threads.size(); ++i) {
			threads[i].join();
		}
	} else {
		std::ostringstream oss;
		for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
			oss << i;

			system(("bash " + oss.str() + Ifeffer::IFEFFIT_SCRIPT).c_str());

			oss.clear();
			oss.str("");
		}
	}
	
}

void Ifeffer::staticEntry(const char* command) {
	system(command);
}

double Ifeffer::run(std::vector<Point> atomic_coordinates, bool threaded) {

	this->updateFEFFFiles(atomic_coordinates);
	this->processIFEFFIT(threaded);

	return this->calculateRMSD();
}

bool Ifeffer::prepare(std::vector<Point> atomic_coordinates, std::string target_exafs_filename, std::string calculated_exafs_filename, double x_min, double x_max) {

	std::string feff_setup_steps[] = {
		"feff6 > /dev/null",
		"for i in `ls feff00*`",
		"do",
		"\techo \"$i\" >> process_files.dat",
		"done",
		"cd .."
	};

	std::string ifeffit_header[] = {
		"read_data(file=chi.chi3,type=chi,group=data)",
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
		"ifeffit -q process.iff"
	};

	std::string setup_feff_script("feff_setup.sh");
	std::ofstream feff_setup(setup_feff_script.c_str());

	// Create folders for each feff/ifeffit run
	
	int size = (int)this->target_indexes.size();
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

		std::ofstream ifeffit_script_file((oss.str() + Ifeffer::IFEFFIT_SCRIPT).c_str());
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
	this->updateFEFFFiles(atomic_coordinates);

	system(("bash " + setup_feff_script).c_str());
	system(("rm " + setup_feff_script).c_str());

	for (int i = 0; i < size; ++i) {
		std::ostringstream oss;
		oss << i;

		std::ifstream process_file((oss.str() + "/process_files.dat").c_str());
		if (!process_file.is_open()) {
			return false;
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
			return false;
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
	this->calculated_exafs_filename = calculated_exafs_filename;
	this->x_min = x_min;
	this->x_max = x_max;

	this->clean_script();

	return true;
}

double Ifeffer::calculateRMSD() {

	// Read in target_indexes output files.
	// Average all files.
	std::vector< std::vector< std::pair<double, double> > > calculated_exafs;
	std::ostringstream oss;

	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		std::vector< std::pair<double, double> > exafs_data;
		oss << i;
		std::ifstream exafs_file(oss.str() + "/" + this->calculated_exafs_filename);

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

	std::vector< std::pair<double, double> > averaged_data;

	// Preload with 0's
	for (int i = 0; i < (int)calculated_exafs[0].size(); ++i) {
		averaged_data.push_back( std::make_pair(calculated_exafs[0][i].first, 0) );
	}

	for (int i = 0; i < (int)calculated_exafs.size(); ++i) {
		for (int j = 0; j < (int)calculated_exafs[i].size(); ++j) {
			averaged_data[j].second += calculated_exafs[i][j].second;
		}
	}

	for (int i = 0; i < (int)averaged_data.size(); ++i) {
		averaged_data[i].second /= calculated_exafs.size();
	}

	std::vector< std::pair<double, double> > calc_data = averaged_data;
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

void Ifeffer::clean() {
	int size = (int)this->target_indexes.size();
	std::ostringstream oss;
	for (int i = 0; i < size; ++i) {
		oss << i;
		system(("rm -rf " + oss.str()).c_str());
		system(("rm -rf " + oss.str() + Ifeffer::IFEFFIT_SCRIPT).c_str());
		oss.clear();
		oss.str("");
	}
}

void Ifeffer::clean_script() {

	std::ofstream cleaner("clean.sh");
	for (int i = 0; i < (int)this->target_indexes.size(); ++i) {
		cleaner << "rm -rf " << i << std::endl;
		cleaner << "rm -rf " << i << Ifeffer::IFEFFIT_SCRIPT << std::endl;
	}
}