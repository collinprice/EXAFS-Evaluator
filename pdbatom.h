#pragma once
#include <string>
#include <string.h>

class PDBAtom {

public:
	
	std::string atomic_symbol;
	double x;
	double y;
	double z;

	PDBAtom(std::string atomic_symbol, int index, double x, double y, double z);
	int getAtomicNumber();
	int getIndex();

	static std::string atomicNumberToSymbol(int atomic_number);
	static int atomicSymbolToNumber(std::string atomic_symbol);

private:

	int index;
	const static int periodic_table_size = 103;
	const static char* periodic_table[];

};