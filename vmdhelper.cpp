#include "vmdhelper.h"

#include <iostream>
#include <fstream>

const std::string VMDHelper::VMD_SCRIPT = "energy_script.vmd";
const std::string VMDHelper::ENERGY_OUTPUT = "pdb_energy";

VMDHelper::VMDHelper(std::string pdb_file, std::string amber_topology_file, std::string namd2_path, std::string vmd_path) {

	std::ofstream vmd_file(VMDHelper::VMD_SCRIPT);
	vmd_file << "package require namdenergy" << std::endl;
	vmd_file << "mol new " << amber_topology_file << " type parm7" << std::endl;
	vmd_file << "mol addfile " << pdb_file << std::endl;
	vmd_file << "set sel [atomselect top \"occupancy 1.0\"]" << std::endl;
	vmd_file << "namdenergy -all -sel $sel -exe " << namd2_path << " -ofile " << VMDHelper::ENERGY_OUTPUT << std::endl;
	vmd_file << "quit" << std::endl;
	vmd_file.close();

	this->vmd_path = vmd_path;
}

float VMDHelper::calculateEnergy() {

	system((this->vmd_path + " -dispdev text -e " + VMDHelper::VMD_SCRIPT + " > /dev/null").c_str());

	std::string line;
	std::ifstream energy_file(VMDHelper::ENERGY_OUTPUT);
	std::getline(energy_file, line); // Skip header line.

	// Kind of a hack but I know "Total" is in the 11th column.
	for (int i = 0; i < 11; ++i) {
		energy_file >> line;
	}

	return atof(line.c_str());
}