#include "chromosome.h"

Chromosome::Chromosome() {

}

Chromosome::Chromosome(std::vector<PDBAtom> points) {

	for (std::vector<PDBAtom>::iterator point = points.begin(); point != points.end(); ++point) {
		
		this->push_back(*point);
	}
}