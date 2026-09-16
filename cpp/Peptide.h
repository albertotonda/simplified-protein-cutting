// class modeling a peptide/protein in the mix; it's basically an extension of std::string, with a few extra methods
// to manage disulfide bonds and other peculiarities of the proteins

#ifndef __PEPTIDE__
#define __PEPTIDE__

#include <iostream>
#include <string>
#include <vector>

class Peptide
{
public:
	// vector with all the disulfide bonds
	std::vector<unsigned int> disulfideBonds;
	// internal string with the peptide
	std::string peptide;

public:
	// constructors/destructor
	Peptide();
	Peptide(std::string peptide);
	~Peptide();
	
	// obtain the substring, but also report information on the bonds inside the two substrings
	Peptide substr(unsigned int start, unsigned int length);
	Peptide substr(unsigned int pos);
	
	// check: is the current position somehow linked to a disulfide bond? 
	bool isCuttable(unsigned int p1, unsigned int p2);
	
	// add a disulfide bond, but first check whether it's valid
	void addDisulfideBond(unsigned int pos);
	
	// length
	inline unsigned int length()
	{
		return this->peptide.length();
	}
	
	// comparison
	int compare(Peptide p)
	{
		return this->peptide.compare( p.peptide );
	}
	
	// access single character
	char& operator[] (size_t p)
	{
		return this->peptide[p];
	}
	
	// compare two peptides
	// NOTE: both "const" here are required, not stylistic -- std::unique() (used on a
	// vector<Peptide> in EndoproteaseModel::run()) compares through const references
	// internally, per the standard's requirements for EqualityComparable. libstdc++
	// (MinGW, Linux) tolerated the previous non-const signature; libc++ (macOS/Clang)
	// correctly rejects it, refusing to compile at all.
	bool operator==(const Peptide& p) const
	{
		return this->peptide == p.peptide;
	}
	
	// used for cout and stuff
	friend std::ostream& operator<<(std::ostream& os, const Peptide& p) 
	{
		os << p.peptide;
		return os;
	}

}; 

#endif
