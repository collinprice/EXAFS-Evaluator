#include "vmdhelper.h"

#include <iostream>
#include <fstream>

const std::string VMDHelper::VMD_RUNNABLE = "vmd.sh";
const std::string VMDHelper::VMD_SCRIPT = "energy_script.vmd";
const std::string VMDHelper::ENERGY_OUTPUT = "pdb_energy";

VMDHelper::VMDHelper(std::string directory, std::string pdb_file, std::string amber_topology_file, std::string namd2_path, std::string vmd_path) {

	// Create script to run vmd.
	std::ofstream vmd_runnable(directory + "/" + VMDHelper::VMD_RUNNABLE);
	vmd_runnable << "cd " << directory << std::endl;
	vmd_runnable << vmd_path << " -dispdev text -e " << VMDHelper::VMD_SCRIPT << " > /dev/null" << std::endl;
	vmd_runnable.close();

	// Create vmd script.
	std::ofstream vmd_file(directory + "/" + VMDHelper::VMD_SCRIPT);
	vmd_file << "package require namdenergy" << std::endl;
	vmd_file << "mol new " << amber_topology_file << " type parm7" << std::endl;
	vmd_file << "mol addfile " << pdb_file << std::endl;
	vmd_file << "set sel [atomselect top \"occupancy 1.0\"]" << std::endl;
	vmd_file << "namdenergy -all -sel $sel -exe " << namd2_path << " -ofile " << VMDHelper::ENERGY_OUTPUT << std::endl;
	vmd_file << "quit" << std::endl;
	vmd_file.close();

	// Copy amber topology file to directory
	system(("cp " + amber_topology_file + " " + directory + "/").c_str());

	this->vmd_path = vmd_path;
	this->directory = directory;
}

double VMDHelper::calculateEnergy() {

	system(("bash " + this->directory + "/" + VMD_RUNNABLE).c_str());

	std::string line;
	std::ifstream energy_file(this->directory + "/" + VMDHelper::ENERGY_OUTPUT);
	std::getline(energy_file, line); // Skip header line.

	// Kind of a hack but I know "Total" is in the 11th column.
	for (int i = 0; i < 11; ++i) {
		energy_file >> line;
	}

	// Clean up
	system(("rm " + directory + "/" + VMDHelper::ENERGY_OUTPUT).c_str());

	return atof(line.c_str());
}