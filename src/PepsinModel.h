// Class representing the pepsin model
// by Alberto Tonda, 2014 <alberto.tonda@gmail.com>

#ifndef __PEPSINMODEL__
#define __PEPSINMODEL__

#include <map>
#include <string>
#include <vector>

// forward declaration
class Peptide;

class PepsinModel
{

public :
	// attributes
	// probability of cutting, given a certain point
	std::map< std::string, double > cutProbability;
	// alterations to the probabilities
	std::map< std::string, std::map<std::string, double> > alterations;
	// alterations for terminals (stored separately)
	std::map<int, double> terminalAlterations;
	// initial proteins
	std::vector<Peptide> originalProteins;
	// statistics on the proteins
	std::map<std::string, std::map<unsigned int, unsigned int> > statistics;
	// pepsin history
	std::vector<double> pepsinHistory;
	// time2 history TODO maybe it could be stored more efficiently
	std::vector<unsigned int> time2History;
	
	// parameters of the model
	double currentPepsin;
	unsigned int overallLength; // length of all proteins
	unsigned int maxAttempts;
	unsigned int maxAttemptsPerTime;
	double maxDH; // maximum degree of hydrolysis
	unsigned int maxTime;
	bool pepsinAlwaysDying;
	double pepsinDyingRatio;
	unsigned int randomSeed;
	unsigned int t;

	// constructor / destructor
	PepsinModel();
	~PepsinModel();

	// methods
	//double computeCutProbability( std::string key, std::string* protein, unsigned int position );
	double computeCutProbability( std::string key, Peptide* protein, unsigned int position );
	void run();
	void step();

	int readJson( std::string fileName );
	int writeLog( std::string fileName );

};

#endif
