#pragma once
#include "pdbatom.h"

#include <vector>
#include <string>
#include <map>

class IFEFFITHelper {

public:
	IFEFFITHelper(std::string folder_name, std::vector<PDBAtom> atoms, std::string target_atom, std::string target_exafs_filename, double x_min, double x_max, std::string feff_path, std::string ifeffit_path);
	~IFEFFITHelper();
	double run(std::vector<PDBAtom> updated_atoms, bool threaded);
	double run(std::string atoms_file, bool threaded);
	std::vector< std::pair<double, double> > getEXAFSData();
	std::vector< std::pair<double, double> > getTargetEXAFS();
	void clean();
	
private:

	static const std::string CALCULATED_EXAFS_FILENAME;
	static const std::string IFEFFIT_SCRIPT;
	static const std::string FEFF_SCRIPT;
	static const std::string CLEAN_SCRIPT;

	std::string target_exafs_filename;
	std::string folder_name;
	std::vector<int> target_indexes;
	std::map<int, int> unique_atoms;
	std::vector<std::string> cached_header;
	std::vector<std::string> cached_feff_filenames;
	std::vector< std::pair<double, double> > cached_target_exafs;
	std::vector< std::pair<double, double> > averaged_calculated_data;
	double x_min;
	double x_max;

	void readTargetEXAFS(std::string filename);
	bool updateProcessFiles();
	bool generateFEFFFile(std::vector<PDBAtom> atomic_coordinates, int index, std::string filename);
	bool updateFEFFFiles(std::vector<PDBAtom> atomic_coordinates);
	void processIFEFFIT(bool threaded);
	static void staticEntry(const char* command);
	double calculateRMSD();

	double failedAttempt();

	bool canPerformIFEFFITCalculations();
	void removeAllCalculatedEXAFSFiles();
	void generateCleanCalculatedEXAFSFilesScript();

};