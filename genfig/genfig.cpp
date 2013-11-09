#include "genfig.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <algorithm> 
#include <functional> 
#include <cctype>
#include <locale>

Genfig::Genfig(std::string filename) {
	init(filename, ':');
}

Genfig::Genfig(std::string filename, char delim) {
	init(filename, delim);
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
        return ltrim(rtrim(s));
}

void Genfig::init(std::string filename, char delim) {

	this->delimeter = delim;
	this->list_delim = ',';

	std::ifstream ifile(filename.c_str());
	if (ifile.good()) {
		std::string line;
		while(std::getline(ifile, line)) {
			trim(line);
			if (!line.empty() && line[0] != '#' && line.find(delim) != std::string::npos) {
				std::stringstream ss(line);
				std::string key;
				std::string value;

				std::getline(ss, key, delim);
				std::getline(ss, value);

				trim(key);
				trim(value);

				this->config[key] = value;
			}
		}
	}
	ifile.close();
}

std::string Genfig::getString(std::string key) {
	return this->config[key];
}

int Genfig::getInt(std::string key) {
	return atoi(this->config[key].c_str());
}

double Genfig::getDouble(std::string key) {
	return atof(this->config[key].c_str());
}

bool Genfig::getBool(std::string key) {
	return this->config[key].compare("true") == 0;
}

std::vector<std::string> Genfig::getStringList(std::string key) {

	std::vector<std::string> stringList;
	std::stringstream ss(this->getString(key));
	for (std::string token; std::getline(ss, token, this->list_delim); ) {
		trim(token);
		if (token.length() > 0) {
			stringList.push_back(token);
		}
	}

	return stringList;
}

std::vector<int>  Genfig::getIntList(std::string key) {
	return this->convertStringListToInts(this->getStringList(key));
}

std::vector<double> Genfig::getDoubleList(std::string key) {
	return this->convertStringListToDoubles(this->getStringList(key));
}

std::vector<bool> Genfig::getBoolList(std::string key) {
	return this->convertStringListToBooleans(this->getStringList(key));
}

bool Genfig::hasKey(std::string key) {
	return this->config.count(key) == 1;
}

void Genfig::writeToFile(std::string filename) {

	std::ofstream output(filename.c_str());
	for (std::map<std::string, std::string>::iterator i = this->config.begin(); i != this->config.end(); ++i) {
		
		output << i->first << ": " << i->second << std::endl;
	}
	output.close();

}

std::vector<int>  Genfig::convertStringListToInts(std::vector<std::string> strings) {

	std::vector<int> ints;
	for (std::vector<std::string>::iterator string = strings.begin(); string != strings.end(); ++string) {
		ints.push_back(atoi((*string).c_str()));
	}
	return ints;
}
std::vector<double>  Genfig::convertStringListToDoubles(std::vector<std::string> strings) {

	std::vector<double> doubles;
	for (std::vector<std::string>::iterator string = strings.begin(); string != strings.end(); ++string) {
		doubles.push_back(atof((*string).c_str()));
	}
	return doubles;
}
std::vector<bool>  Genfig::convertStringListToBooleans(std::vector<std::string> strings) {

	std::vector<bool> bools;
	for (std::vector<std::string>::iterator string = strings.begin(); string != strings.end(); ++string) {
		bools.push_back( !(*string).compare("true") );
	}
	return bools;
}