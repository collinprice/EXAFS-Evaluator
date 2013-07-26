#include "molecule.h"

#include <vector>
#include <map>
#include <utility>
#include <string>

class Ifeffer {

	private:
		// Atom* target_atom;
		std::vector<Atom> atoms;
		std::vector<int> target_indexes;
		std::map<int, int> unique_atoms;

		std::vector<std::string> cached_header;
		std::vector<std::string> cached_feff_filenames;
		std::vector< std::pair<double, double> > cached_target_exafs;
		double x_min;
		double x_max;

		std::string calculated_exafs_filename;
		static const std::string IFEFFIT_SCRIPT;

		void readTargetEXAFS(std::string filename);

		bool generateFEFFFile(std::vector<Point> atomic_coordinates, int index, std::string filename);
		bool updateFEFFFiles(std::vector<Point> atomic_coordinates);
		void processIFEFFIT(bool threaded);
		static void staticEntry(const char* command);

		double calculateRMSD();

		void clean_script();

	public:

		Ifeffer();
		Ifeffer(std::vector<Atom> atoms, Atom target_atom);

		bool prepare(std::vector<Point> atomic_coordinates, std::string target_exafs_filename, std::string calculated_exafs_filename, double x_min, double y_min);
		double run(std::vector<Point> atomic_coordinates, bool threaded);

		void clean();

};