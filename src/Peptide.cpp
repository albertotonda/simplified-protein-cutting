// methods for the Peptide class
#include "Peptide.h"
#include <Log.h>

using namespace std;

// constructors
Peptide::Peptide() :
peptide("")
{}

Peptide::Peptide(string peptide) :
peptide(peptide)
{}

// the destructor is basically the same, but also deletes the associated vectors
Peptide::~Peptide()
{
	this->peptide.clear();
	this->disulfideBonds.clear();
}

// this method is the basic substr method, but it also reports information on the bonds 
Peptide Peptide::substr(unsigned int start, unsigned int length )
{
	Peptide result;

	// get the substr
	result.peptide = this->peptide.substr(start, length);
	
	// also report the bond information onto the new peptide, normalized on the new length
	unsigned int end = start + length;
	for(unsigned int i = 0; i < this->disulfideBonds.size(); i++)
	{
		// if the disulfide bond is inside the new substring
		if( this->disulfideBonds[i] >= start && this->disulfideBonds[i] < end )
		{
			// recreate the disulfide bond in the product, with the correct offset
			result.disulfideBonds.push_back( this->disulfideBonds[i] - start );
			LOG_TRACE(	"The old peptide had a disulfide bond in position #" << this->disulfideBonds[i]
					<< " (\"" << this->peptide[ this->disulfideBonds[i] ]
					<< "\"), while the new peptide has a disulfide bond in position #"
					<< (this->disulfideBonds[i] - start)
					<< " (\"" << result.peptide[ this->disulfideBonds[i] - start ] << "\")" );
		}
	}

	return result;
}

// overload: from the current position to the end of the string
Peptide Peptide::substr(unsigned int pos)
{
	Peptide result;
	result.peptide = this->peptide.substr(pos);
	
	for(unsigned int i = 0; i < this->disulfideBonds.size(); i++)
	{
		if( this->disulfideBonds[i] >= pos )
		{
			result.disulfideBonds.push_back( this->disulfideBonds[i] - pos );
			LOG_TRACE(	"The old peptide had a disulfide bond in position #" << this->disulfideBonds[i]
					<< " (\"" << this->peptide[ this->disulfideBonds[i] ]
					<< "\"), while the new peptide has a disulfide bond in position #"
					<< (this->disulfideBonds[i] - pos)
					<< " (\"" << result.peptide[ this->disulfideBonds[i] - pos ] << "\")" );
		}
	}
	
	return result;
}

// check whether a certain position is cuttable
bool Peptide::isCuttable(unsigned int p1, unsigned int p2)
{
	// check all disulfide bonds; basically, you can't cut around it
	for(unsigned int i = 0; i < this->disulfideBonds.size(); i++)
	{
		if( p1 == this->disulfideBonds[i] || p2 == this->disulfideBonds[i] )
			return false;
	}

	return true;
}

void Peptide::addDisulfideBond(unsigned int pos)
{
	if( pos < this->peptide.length() )
	{
		for(unsigned int i = 0; i < this->disulfideBonds.size(); i++) 
		{
			if( this->disulfideBonds[i] == pos)
			{
				LOG_WARN("Warning: bond in position #" << pos << " is already inside the list of disulfide bonds.");
				return;
			}
		}

		this->disulfideBonds.push_back( pos );

		LOG_DEBUG("Added new disulfide bond in position #" << pos << ", corresponding to amino-acid \"" << this->peptide[pos] << "\"");

		return;
	}

	LOG_WARN("Warning: cannot add disulfide bond in position #" << pos << ": non-valid position.");
}
