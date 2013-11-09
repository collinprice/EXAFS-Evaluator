#include <string>
#include <map>
#include <vector>

class Genfig {

	public:
		Genfig(std::string filename);
		Genfig(std::string filename, char delim);
		std::string getString(std::string key);
		int getInt(std::string key);
		double getDouble(std::string key);
		bool getBool(std::string key);

		std::vector<std::string> getStringList(std::string key);
		std::vector<int> getIntList(std::string key);
		std::vector<double> getDoubleList(std::string key);
		std::vector<bool> getBoolList(std::string key);

		bool hasKey(std::string key);
		void writeToFile(std::string filename);

	private:
		std::map<std::string, std::string> config;
		char delimeter;
		char list_delim;

		std::vector<int> convertStringListToInts(std::vector<std::string> strings);
		std::vector<double> convertStringListToDoubles(std::vector<std::string> strings);
		std::vector<bool> convertStringListToBooleans(std::vector<std::string> strings);

		void init(std::string filename, char delim);
		
};