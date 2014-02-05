#pragma once
#include <string>

class VMDHelper {

public:
	VMDHelper(std::string directory, std::string pdb_file, std::string amber_topology_file, std::string namd2_path, std::string vmd_path);

	double calculateEnergy();

private:

	static const std::string VMD_RUNNABLE;
	static const std::string VMD_SCRIPT;
	static const std::string ENERGY_OUTPUT;
	std::string vmd_path;
	std::string directory;
};