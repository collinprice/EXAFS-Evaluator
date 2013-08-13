#include "clustering.h"
#include <iostream>

std::vector< std::vector<double> > Clustering::createTable(std::vector< std::vector<PDBAtom> > atoms) {

	std::vector< std::vector<double> > table;

	for (int i = 0; i < (int)atoms.size(); ++i) {

		std::cout << i << std::endl;

		std::vector<double> row;
		for (int j = 0; j < (int)atoms.size(); ++j) {

			double difference = 0;
			for (int k = 0; k < (int)atoms[0].size(); ++k) {
				
				difference += atoms[i][k].distance(atoms[j][k]);
			}
			row.push_back(difference);

			if (i == j) break;
		}

		table.push_back(row);
	}

	return table;
}