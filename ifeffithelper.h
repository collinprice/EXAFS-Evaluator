#include "pdbatom.h"

#include <vector>
#include <string>
#include <map>

class IFEFFITHelper {

public:
	IFEFFITHelper(std::vector<PDBAtom> atoms, std::string target_atom, std::string target_exafs_filename, double x_min, double x_max, std::string feff_path, std::string ifeffit_path);
	~IFEFFITHelper();
	double run(std::vector<PDBAtom> updated_atoms, bool threaded);
	std::vector< std::pair<double, double> > getEXAFSData();
	std::vector< std::pair<double, double> > getTargetEXAFS();
	void clean();
	
private:

	static const std::string CALCULATED_EXAFS_FILENAME;
	static const std::string IFEFFIT_SCRIPT;
	static const std::string CLEAN_SCRIPT;

	std::vector<int> target_indexes;
	std::map<int, int> unique_atoms;
	std::vector<std::string> cached_header;
	std::vector<std::string> cached_feff_filenames;
	std::vector< std::pair<double, double> > cached_target_exafs;
	std::vector< std::pair<double, double> > averaged_calculated_data;
	double x_min;
	double x_max;

	void readTargetEXAFS(std::string filename);
	bool generateFEFFFile(std::vector<PDBAtom> atomic_coordinates, int index, std::string filename);
	bool updateFEFFFiles(std::vector<PDBAtom> atomic_coordinates);
	void processIFEFFIT(bool threaded);
	static void staticEntry(const char* command);
	double calculateRMSD();
	void clean_script();

	bool canPerformIFEFFITCalculations();
	void removeAllCalculatedEXAFSFiles();

};