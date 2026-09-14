// Methods for the class representing the pepsin model

// class header
#include <PepsinModel.h>

// other classes used by this model
#include <Peptide.h>

// standard classes
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>

// local classes/libraries
#include <tinyxml.h>

// some defines, used mainly for the XML parsing
#define DEFAULT_CSV "statistics-"
#define DEFAULT_OUTPUT "solution.xml"
#define DEFAULT_REPETITIONS 1

#define MAX_ATTEMPTS 10
#define MAX_FAILED_ATTEMPTS 100
#define MAX_TIME 10000

#define NAME_PEPSIN "pepsin"

#define XML_ALTERATION "alteration"
#define XML_ALTERATIONS "alterations"
#define XML_CONFIGURATION "configuration"
#define XML_CUT "cut"
#define XML_CUTS "cuts"
#define XML_INITIALPEPSIN "initialPepsin"
#define XML_MAXATTEMPTS "maxAttempts"
#define XML_MAXATTEMPTSPERTIME "maxAttemptsPerTime"
#define XML_MAXDH "maxDH"
#define XML_MAXTIME "maxTime"
#define XML_PARAMETERS "parameters"
#define XML_PEPSINALWAYSDYING "pepsinAlwaysDying"
#define XML_PEPSINDYINGRATIO "pepsinDyingRatio"
#define XML_PROTEIN "protein"
#define XML_PROTEINS "proteins"
#define XML_RANDOMSEED "randomSeed"
#define XML_TERMINAL "terminal"

#define XML_ATTRIBUTE_AMINOACID "aminoacid"
#define XML_ATTRIBUTE_LEFT "left"
#define XML_ATTRIBUTE_MULTIPLIER "multiplier"
#define XML_ATTRIBUTE_NAME "name"
#define XML_ATTRIBUTE_RIGHT "right"
#define XML_ATTRIBUTE_POSITION "position"
#define XML_ATTRIBUTE_DISULFIDEBONDS "disulfideBonds"
#define XML_ATTRIBUTE_PROBABILITY "probability"
#define XML_ATTRIBUTE_QUANTITY "quantity"
#define XML_ATTRIBUTE_SEED "randomSeed"
#define XML_ATTRIBUTE_SIDE "side"
#define XML_ATTRIBUTE_VALUE "value"

// debugging purposes
#define DEBUG 1 // set to 1 to activate debugging

// namespace
using namespace std;

// constructor (basically it does nothing at all, besides providing a few default values)
PepsinModel::PepsinModel() :
currentPepsin(1.0),
maxAttemptsPerTime(10),
maxAttempts(1000),
maxDH(100.0),
maxTime(10000),
overallLength(0),
pepsinAlwaysDying(false),
pepsinDyingRatio(1.0),
randomSeed(0),
t(0),
verbose(false)
{}

// destructor
PepsinModel::~PepsinModel()
{
	// clear all the maps and vectors
	this->alterations.clear();
	this->cutProbability.clear();
	this->originalProteins.clear();
	this->pepsinHistory.clear();
	this->statistics.clear();
}

// read the XML file and store the relevant information 
int PepsinModel::readXml( string fileName )
{
	if(verbose) cout << "Loading XML file..." << endl;
	TiXmlDocument doc(fileName.c_str());

	if( !doc.LoadFile() )
	{
		cerr << "Error: cannot open file \"" << fileName << "\", or file is not proper XML. Aborting..." << endl;
		return -1;
	}
	
	// start with the root
	if(verbose) cout << "Parsing XML file..." << endl;
	TiXmlElement* pRoot = doc.FirstChildElement(XML_CONFIGURATION);
	if( pRoot == NULL )
	{
		cerr << "Error: tag \"" << XML_CONFIGURATION << "\" not found in the document. Aborting..." << endl;
		return -1;
	}
	
	// first of all, parse the "parameters"
	TiXmlElement* pParameters = pRoot->FirstChildElement(XML_PARAMETERS);
	
	// parameters are not compulsory (TODO REALLY?)
	
	if( pParameters != NULL )
	{
		// check if there's a random seed
		TiXmlElement* pRandomSeed = pParameters->FirstChildElement(XML_RANDOMSEED);
		if( pRandomSeed != NULL )
		{
			const char* valueString = pRandomSeed->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%u", &this->randomSeed);
		}
		
		// check if there is a maximum time specified
		TiXmlElement* pMaxTime = pParameters->FirstChildElement(XML_MAXTIME);
		if( pMaxTime != NULL )
		{
			const char* valueString = pMaxTime->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%u", &this->maxTime);
		}

		// check if there is a maximum time specified
		TiXmlElement* pMaxDH = pParameters->FirstChildElement(XML_MAXDH);
		if( pMaxDH != NULL )
		{
			const char* valueString = pMaxDH->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%lf", &this->maxDH);
		}

		// check for maximum attempts per time 
		TiXmlElement* pMaxAttemptsPerTime = pParameters->FirstChildElement(XML_MAXATTEMPTSPERTIME);
		if( pMaxAttemptsPerTime != NULL )
		{
			const char* valueString = pMaxAttemptsPerTime->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%u", &this->maxAttemptsPerTime);
		}

		// check for maximum attempts
		TiXmlElement* pMaxAttempts = pParameters->FirstChildElement(XML_MAXATTEMPTS);
		if( pMaxAttempts != NULL )
		{
			const char* valueString = pMaxAttempts->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%u", &this->maxAttempts);
		}

		// check for initial pepsin
		TiXmlElement* pInitialPepsin = pParameters->FirstChildElement(XML_INITIALPEPSIN);
		if( pInitialPepsin != NULL )
		{
			const char* valueString = pInitialPepsin->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%lf", &this->currentPepsin);
		}
		
		// check if pepsin dies at every timestep, or just when it does not cut
		TiXmlElement* pPepsinAlwaysDying = pParameters->FirstChildElement(XML_PEPSINALWAYSDYING);
		if( pPepsinAlwaysDying != NULL )
		{
			const char* valueString = pInitialPepsin->Attribute(XML_ATTRIBUTE_VALUE);
			if( strcmp(valueString, "false") == 0 || strcmp(valueString, "0") == 0 )
				this->pepsinAlwaysDying = false;
			else if( strcmp(valueString, "true") == 0 || strcmp(valueString, "1") == 0 )
				this->pepsinAlwaysDying = true;
			else
				cerr 	<< "Warning: value \"" << valueString 
					<< "\" is not a valid option for tag \"" << XML_PEPSINALWAYSDYING
					<< "\" (valid values are \"true\"/1 and \"false\"/0). The pepsinAlwaysDying will be set to "
					<< this->pepsinAlwaysDying
					<< endl;
		}		

		// check for pepsin dying ratio
		TiXmlElement* pPepsinDyingRatio = pParameters->FirstChildElement(XML_PEPSINDYINGRATIO);
		if( pPepsinDyingRatio != NULL )
		{
			const char* valueString = pPepsinDyingRatio->Attribute(XML_ATTRIBUTE_VALUE);
			sscanf(valueString, "%lf", &this->pepsinDyingRatio);
		}

	}

	if( randomSeed != 0 )
	{
		// initialize random generator with seed
		if(verbose) cout << "Initializing random generator with seed " << randomSeed << endl;
		srand( randomSeed );
	}
	else
	{
		// initialize random generator with t
		time_t timeRandomSeed = time(NULL);
		if(verbose) cout << "Initializing random generator with time (" << timeRandomSeed << ")" << endl;
		srand( timeRandomSeed );
	}
	
	// take the proteins
	// include the possibility of having different proteins in the initial set
	TiXmlElement* pProtein = pRoot->FirstChildElement(XML_PROTEINS);
	if( pProtein == NULL )
	{
		cerr << "Error: tag \"" << XML_PROTEINS << "\" must be specified. Aborting..." << endl;
		return -1;
	}
	
	// go into the tag
	pProtein = pProtein->FirstChildElement(XML_PROTEIN);
	
	while( pProtein != NULL )
	{
		if(verbose) cout	<< "Adding a protein..." << endl;
		string originalProtein = pProtein->GetText();
		
		// ok, here we get two possible attributes: the name of the protein (unused) and the quantity (important!)
		// recently, we added a third attribute, a string with a list of disulfide bonds
		const char* proteinName = pProtein->Attribute(XML_ATTRIBUTE_NAME);
		if(proteinName != NULL && verbose) cout	<< "Adding protein \"" << proteinName << "\"..." << endl;

		// default quantity for each protein is 1
		unsigned int proteinQuantity = 1;
		const char* proteinQuantityString = pProtein->Attribute(XML_ATTRIBUTE_QUANTITY);
		if( proteinQuantityString != NULL )
		{
			sscanf(proteinQuantityString, "%u", &proteinQuantity);
		}
		
		// remove all whitespaces and other stuff from the string
		originalProtein.erase( std::remove_if( originalProtein.begin(), originalProtein.end(), ::isspace ), originalProtein.end() );
		
		if(verbose) cout 	<< "Adding protein with composition \"" << originalProtein 
					<< "\", " << originalProtein.length() 
					<< " amminoacids long."  << endl;
		
		// create the Peptide with the original protein
		Peptide originalProteinPeptide(originalProtein);

		// add disulfide bonds (if they're there)
		const char* disulfideBonds = pProtein->Attribute(XML_ATTRIBUTE_DISULFIDEBONDS);
		if(disulfideBonds != NULL)
		{
			if(verbose) cout 	<< "Found list of disulfide bonds: \"" << disulfideBonds << "\"" << endl;

			// parse the string to get a list of unsigned int to add to the Peptide class
			stringstream stream( disulfideBonds );
			do
			{
				unsigned int bond;
				stream >> bond;
				
				if(verbose) cout 	<< "Adding bond in position #" << bond
							<< " to the protein (modified to #" 
							<< (bond-1) << " for C++ internal array indexing)"
							<< endl;

				originalProteinPeptide.addDisulfideBond( bond-1 );
			}while(stream);
			
		}
		
		for(unsigned int i = 0; i < proteinQuantity; i++)
		{
			if(verbose) cout << "\tAdding copy #" << (i+1) << "..." << endl;
			// add the protein to the initial set
			this->originalProteins.push_back( originalProteinPeptide );
			
			// increase the overall length of the original proteins
			this->overallLength += originalProtein.length();
		}
		
		// on to the next protein
		pProtein = pProtein->NextSiblingElement();
	}
	
	// check if there are proteins at all
	if( originalProteins.size() == 0 )
	{
		cerr << "Error: no tags \"" << XML_PROTEIN << "\" found. Aborting..." << endl;
		return -1;
	}

	// then, take all the cuts
	if(verbose) cout << "Reading the cuts..." << endl;
	TiXmlElement* pCut = pRoot->FirstChildElement(XML_CUTS);
	
	if( pCut == NULL )
	{
		cerr << "Error: tag \"" << XML_CUTS << "\" must be specified. Aborting..." << endl;
		return -1;
	}
	
	// go into the tag
	pCut = pCut->FirstChildElement(XML_CUT);
	
	while( pCut != NULL )
	{
		const char* left;
		const char* right;
		double probability = 0.0;

		// get information about the cut
		left = pCut->Attribute(XML_ATTRIBUTE_LEFT);
		right = pCut->Attribute(XML_ATTRIBUTE_RIGHT);
		pCut->QueryDoubleAttribute(XML_ATTRIBUTE_PROBABILITY, &probability);

		// put everything inside the map, if probability > 0
		if( probability > 0.0 )
		{
			string key = left;
			key += "\t";
			key += right;
			this->cutProbability[ key ] = probability;
			
			if(verbose) cout << "Cut \"" << key << "\" with probability=" << probability << " added" << endl;
		}

		// next element
		pCut = pCut->NextSiblingElement();
	}
	
	// check if the map is empty, it could be a problem
	if( cutProbability.size() == 0 )
	{
		cerr << "Warning: there are no cuts specified, so probably nothing will happen during this simulation..." << endl;
	}
	else
	{
		if(verbose) cout << "There are " << cutProbability.size() << " cuts specified" << endl;
	}

	// proceed with managing the alterations
	TiXmlElement* pAlteration = pRoot->FirstChildElement(XML_ALTERATIONS);
	
	// alterations are optional, so if they're not there it's not a big deal
	if( pAlteration != NULL )
	{
		if(verbose) cout << "Now reading alterations..." << endl;
		pAlteration = pAlteration->FirstChildElement(XML_ALTERATION);
		while( pAlteration != NULL )
		{
			// so, the idea is to put the alterations in a hash map, that initially contains the
			// unaltered probabilities, then compute the average + std influence on each column, 
			// and then evaluate how much each element is outside some standard deviations from
			// the mean; it could actually be a double hash map
			
			// read the attributes
			string aminoacid, side, position, key;
			double probability;
			
			aminoacid = pAlteration->Attribute(XML_ATTRIBUTE_AMINOACID);
			side = pAlteration->Attribute(XML_ATTRIBUTE_SIDE);
			position = pAlteration->Attribute(XML_ATTRIBUTE_POSITION);
			pAlteration->QueryDoubleAttribute(XML_ATTRIBUTE_PROBABILITY, &probability);
			
			if( side.compare("left") == 0 ) 
				key = "-";
			else 
				key = "+";
			key += position;
			
			// store the alteration inside the map
			this->alterations[ key ][ aminoacid ] = probability;
			if(verbose) cout 	<< "- alteration[ " << key << " ][ " << aminoacid << " ] = " 
						<< this->alterations[key][aminoacid] << endl;

			// onto the next alteration
			pAlteration = pAlteration->NextSiblingElement(XML_ALTERATION);
		}
		
		// we added a second type of alteration, called "terminal", to take into account
		// the alteration of the probabilities of cutting, if you are at either end of a chain
		if(verbose) cout << "Now reading terminal alterations..." << endl;
		TiXmlElement* pTerminal = pRoot->FirstChildElement(XML_ALTERATIONS)->FirstChildElement(XML_TERMINAL);
		while( pTerminal != NULL )
		{
			// read the attributes
			string side;
			int key;
			double multiplier;
			
			// TODO error control? raise exceptions?	
			side = pTerminal->Attribute(XML_ATTRIBUTE_SIDE);
			pTerminal->QueryIntAttribute(XML_ATTRIBUTE_POSITION, &key);
			pTerminal->QueryDoubleAttribute(XML_ATTRIBUTE_MULTIPLIER, &multiplier);
			
			if( side.compare("left") == 0 ) 
				key = -1 * key;
			
			// store this alteration inside the map
			this->terminalAlterations[ key ] = multiplier;
			if(verbose) cout	<< "- terminal alteration[ " << key << " ] = "
						<< this->terminalAlterations[key] << endl;
			
			// onto to the next alteration
			pTerminal = pTerminal->NextSiblingElement(XML_TERMINAL);
		}
		
		
		// some statistics on alterations here
		// iterate over the keys to get the average influence
		for( map< string, map<string, double> >::iterator it = this->alterations.begin(); it != this->alterations.end(); it++)
		{
			double average = 0.0;
			for( map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++ )
			{
				average += it2->second;
			}
			average /= it->second.size();
			
			double std = 0.0;
			for( map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++ )
			{
				double value = it2->second - average;
				if( value < 0 ) value = -1.0 * value;

				std += value;
			}
			std /= it->second.size();
			
			if(verbose)	cout 	<< "For alterations on column \"" << it->first 
						<< "\", avg=" << average 
						<< ", std=" << std 
						<< ", sigma=" << sqrt(std)
						<< endl;
			
			// remake the hash map, inserting only the values that will actually modify the probabilities;
			// the alterations might be 
			for( map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++ )
			{
				double value = it2->second - average;
				if( value < 0 ) value = -1.0 * value;

				if( value > std )
				{
					if(verbose) cout 	<< "- alteration[ " << it->first << " ][ " 
								<< it2->first << " ]=" << it2->second 
								<< ", difference from average is " << (value / std ) << " variances!"
								<< endl;
				}
				
				if( value > sqrt(std) )
				{
					if(verbose) cout 	<< "- alteration[ " << it->first << " ][ " 
								<< it2->first << " ]=" << it2->second 
								<< ", difference from average is " << (value / sqrt(std) ) << " standard deviations!"
								<< endl;
				}
				
				it2->second /= average;
			}
		}
		
		
		// this is just debugging / pre-processing, and it should probably be removed; but it is interesting
		// to understand how the cut probabilities are changed by the alterations: so we iterate over the
		// proteins and check the difference with the proposed approach
		if(verbose)
		for(unsigned int p = 0; p < 1 /* originalProteins.size() */; p++)
		{
			cout << "Analyzing the influence of alterations for protein \"" << originalProteins[p] << "\"..." << endl;
			for(int i = 0; i < originalProteins[p].length(); i++)
			{
				// check the cut probability for a couple of amminoacids
				string key = ""; 
				key += originalProteins[p][i];
				key += "\t";

				if( i + 1 < originalProteins[p].length() )
					key += originalProteins[p][i+1];
				
				// verify if there is a corresponding entry in the hash table
				map<string,double>::iterator it = cutProbability.find( key );
				if( it != cutProbability.end() )
				{
					double probability = it->second;
					cout << "- \"" << key << "\": original cut probability = " << probability << endl;
					
					probability = this->computeCutProbability( key, &originalProteins[p], i );
					cout << "- \"" << key << "\": altered cut probability = " << probability << endl;

					cout << endl; // for readability
				}
				
			}
		}
		
		
	} // end if there are alterations
	
	
	// everything all right!
	return 0;
}

// write the history to a file
int PepsinModel::writeLog( string fileName )
{
	if(this->verbose) cout << endl << "Pre-processing statistics..." << endl;

	// trying to put some buffering, in order not to block the file for too long
	stringstream outStream;

	// iterate over each protein to write the header
	outStream << "\"time\",\"time2\",\"" << NAME_PEPSIN << "\"";
	for(map<string, map<unsigned int, unsigned int> >::iterator 	it = statistics.begin(); 
									it != statistics.end(); 
									it++)
	{
		outStream << ",\"" << it->first << "\"";
	}
	outStream << endl;

	// create a temporary structure with all peptides, used later to store the last quantity
	// appearing in the table; many lines are the same, and will be copied
	map<string, unsigned int> lastQuantity;
	
	for(map<string, map<unsigned int, unsigned int> >::iterator 	it = statistics.begin(); 
									it != statistics.end(); 
									it++)
	{
		if( it->second.find( 0 ) != it->second.end() )
			lastQuantity[ it->first ] = it->second.find( 0 )->second;
		else
			lastQuantity[ it->first ] = 0;
	}
	
	// then, iterate over t; if there is no quantity of protein for that time, put 0
	// statistics are stored only "period" iterations, to avoid HUGE log files
	// TODO
	// - is it possible to change the increments so that localt += period?
	// - period readable by configuration file
	unsigned int period = 10; // when localt % period is used, period = 100
	int lastTime2 = -1;
	for(unsigned int localt = 0; localt < t; localt++)
	{
		//if( localt % period == 0 ) outStream 	<< localt << "," 
		if( this->time2History[localt] % period == 0 && this->time2History[localt] != lastTime2 ) 
		{
			outStream 	<< localt << "," 
					<< this->time2History[localt] << ","
					<< this->pepsinHistory[localt];
		}

		for(map<string, map<unsigned int, unsigned int> >::iterator 	it = statistics.begin(); 
										it != statistics.end(); 
										it++)
		{
			if( this->time2History[localt] % period == 0 && this->time2History[localt] != lastTime2 ) 
				outStream << ",";

			// if there is an occurrence in the map for that protein
			if( it->second.find(localt) != it->second.end() )
			{
				if( this->time2History[localt] % period == 0 && this->time2History[localt] != lastTime2 ) 
					outStream << it->second.find(localt)->second;
				
				// update the "last quantity" map
				lastQuantity[ it->first ] = it->second.find(localt)->second;
			}
			else
			{
				// otherwise, the quantity is the last one stored in the map 
				if( this->time2History[localt] % period == 0 && this->time2History[localt] != lastTime2 ) 
					outStream << lastQuantity[ it->first ];
			}
		}

		if( this->time2History[localt] % period == 0 && this->time2History[localt] != lastTime2 ) 
		{
			// end the current line in the buffer
			outStream << endl;
			
			// also, store the new value for t2
			lastTime2 = this->time2History[localt];
		}
	}

	if( this->verbose ) cout 	<< "Writing statistics to CSV file \"" 
					<< fileName << "\"..." 
					<< endl;
	// open file
	ofstream csvOut( fileName.c_str() );
	if( !csvOut.is_open() )
	{
		cerr << "Error: cannot write on file \"" << fileName << "\". Aborting..." << endl;
		return -1; 
	}
	
	csvOut << outStream.str();

	csvOut.close();

	// everything all right!
	return 0;
}

// after initialization, here is the "true" run of the model!
void PepsinModel::run()
{
	cout 	<< "Starting the simulation, with maxTime=" << maxTime

		<< ", maxAttemptsPerTime=" << maxAttemptsPerTime
		<< ", maxAttempts=" << maxAttempts
		<< ", maxDH=" << maxDH
		<< ", initialPepsin=" << currentPepsin
		<< ", pepsinDyingRatio=" << pepsinDyingRatio
		<< endl;

	// initialize some values
	unsigned int time2 = 0; // this variable is increased only when there is a cut
	unsigned int attempts = 0;
	double dh = 0; // current degree of hydrolysis 
	double dhTotalSites = 0;
	//vector<string> currentProteins = this->originalProteins;
	vector<Peptide> currentProteins = this->originalProteins;

	// store the statistics for the current protein; besides t = 0, the following stats
	// are derived by difference with the previous ones; also, create map with last known quantity
	map<string, unsigned int> lastQuantity;
	
	//vector<string> uniqueProteins = currentProteins;
	vector<Peptide> uniqueProteins = currentProteins;
	uniqueProteins.erase( unique( uniqueProteins.begin(), uniqueProteins.end() ), uniqueProteins.end() );

	for(unsigned int i = 0; i < uniqueProteins.size(); i++)
	{
		unsigned int occurrences = 0;
		for(unsigned int j = 0; j < currentProteins.size(); j++)
		{
			if( currentProteins[j].compare( uniqueProteins[i] ) == 0 )
				occurrences++;
		}

		statistics[uniqueProteins[i].peptide][t] = occurrences;
		lastQuantity[uniqueProteins[i].peptide] = occurrences;

		// in order to compute the current degree of hydrolysis, we also need to take into account
		// the total number of sites where the pepsin can cut the proteins
		dhTotalSites += (uniqueProteins[i].length() - 1) * occurrences;
	}
	
	// some debugging
	if( verbose ) cout << "The total number of sites that the pepsin can cut is " << dhTotalSites << endl;
	
	// also, take note of pepsin quantity and time2
	this->pepsinHistory.push_back( currentPepsin );
	this->time2History.push_back( time2 );

	// start!
	while( this->t < this->maxTime && attempts < this->maxAttempts && (time2/(double)dhTotalSites) < this->maxDH )
	{
		cout 	<< endl 
			<< "Time #" << t 
			<< " (Time2 #" << time2 
			<< ", DH=" << (time2/(double)dhTotalSites) 
			<< ", maxDH=" << this->maxDH
			<< ")" <<  endl;
		
		// some variables need to be defined here
		bool cutPerformed = false;
		//string parentProtein, product1, product2;
		Peptide parentProtein, product1, product2;

		// reset the attempts done per instant of t
		unsigned int attemptsPerTime = 0;

		// before doing all the rest, just verify whether the pepsin acts
		double randomPepsinActivation = (double) rand() / RAND_MAX;
		
		if( randomPepsinActivation < currentPepsin)
		{
			// find a random point in a random protein
			// TODO this part could be sped up by creating a map< unsigned int, pair<unsigned int, unsigned int> > 
			//	that keeps track of the protein and relative position for an absolute position;
			//	so, if I have 200 peptides with 20 aminoacids each, I have a total of 4000 positions;
			//	the map can tell me that position 201 is actually position 0 on protein 1

			// ok, *probably* the probability of a random protein to be chosen is 
			// dependent on its size; so, we should probably just choose a random point in the starting protein and
			// go over all the proteins until you find the corresponding point
			unsigned int startingPosition = rand() % overallLength;
			
			bool positionFound = false;
			unsigned int proteinsIndex = 0;
			while( positionFound != true && proteinsIndex < currentProteins.size() )
			{
				if( startingPosition >= currentProteins[proteinsIndex].length() )
				{
					startingPosition -= currentProteins[proteinsIndex].length();
					proteinsIndex++;
				}
				else
				{
					positionFound = true;
				}
			}
			unsigned int startingProtein = proteinsIndex;
			
			if(verbose) cout 	<< "Starting position is " << startingPosition 
						<< " in protein #" << startingProtein 
						<< " (\"" << currentProteins[startingProtein] << "\")"
						<< endl;

			// TODO check for errors, but it should be unlikely...
			
			// now, it's enough to consider the two characters
			// to generalize, you could just check which are the numbers of characters in the cuts
			unsigned int currentPosition = startingPosition;

			while( cutPerformed == false && attemptsPerTime < this->maxAttemptsPerTime )
			{
				string key = ""; 
				key += currentProteins[proteinsIndex][currentPosition];
				key += "\t";

				if( currentPosition + 1 < currentProteins[proteinsIndex].length() )
					key += currentProteins[proteinsIndex][currentPosition+1];
				
				if(verbose) cout 	<< "- Analyzing protein[" << proteinsIndex << "][" << currentPosition 
							<< "]=\"" << key << "\"..." << endl;

				// obtain the probability of cutting in that point
				double probability = this->computeCutProbability( key, &currentProteins[proteinsIndex], currentPosition );

				if( probability > 0.0 )
				{
					// random number between 0 and 1
					double randomNumber = (double) rand() / RAND_MAX;
						
					if( randomNumber < probability )
					{
						// perform the cut!
						if(verbose) cout	<< "Performing cut \"" << key << "\" in position " 
									<< currentPosition << " that had probability " 
									<< cutProbability[key] 
									<< " (modified to " << probability 
									<< ", randomNumber was " << randomNumber << ")"
									<< endl;
						
						// now, the current proteins will change! add the resulting sub-proteins
						// at the end of the vector, and remove the current protein
						product1 = currentProteins[proteinsIndex].substr(0, currentPosition+1);
						product2 = currentProteins[proteinsIndex].substr(currentPosition+1);
						
						if(verbose) cout 	<< "Protein \"" << currentProteins[proteinsIndex] 
									<< "\" has been cut in position " << currentPosition
									<< " and the two resulting products are \"" << product1
									<< "\" and \"" << product2 << "\""
									<< endl;

						// remove parent protein, but keep its structure for the statistics
						parentProtein = currentProteins[proteinsIndex];
						currentProteins.erase( currentProteins.begin() + proteinsIndex );
						
						// add children proteins
						currentProteins.push_back( product1 );
						currentProteins.push_back( product2 );

						// set the flag
						cutPerformed = true;
						
						// also, every time a cut succeeds, reset the number of failed attempts
						attemptsPerTime = 0;
						attempts = 0;
					}
					else
					{
						if(verbose) cout 	<< "Probability was " << probability
									<< ", while randomNumber was " << randomNumber
									<< ": cut not performed."
									<< endl;
					}
					
				}
				else
				{
					if(verbose) cout 	<< "Probability to cut = 0.0; there are no cuts for this key, or there's a disulfide bond near the target link." << endl;
				}
				// if there is no cut available for the two letters, skip to the next position
				currentPosition++;
				
				// if the position is out of bounds
				if( cutPerformed == false && currentPosition >= currentProteins[proteinsIndex].length() )
				{
					proteinsIndex++;
				
					// if we go over the maximum number of proteins
					if( proteinsIndex >= currentProteins.size() )
					{
						proteinsIndex = 0;
						if(verbose) cout 	<< "Going over the list of proteins again: attemptsPerTime=" 
									<< attemptsPerTime << ", attempts=" << attempts	
									<< ", startingPosition=" << startingPosition
									<< ", startingProtein=" << startingProtein
									<< endl;
					}
					

					currentPosition = 0;
					if(verbose) cout 	<< "- Now analyzing protein #" << proteinsIndex 
								<< " (size " << currentProteins[proteinsIndex].length() << ")"
								<< endl;
				}

				// if we came back where we started...
				if( cutPerformed == false )
				{
					attemptsPerTime++;
					attempts++;
					
					// also, reduce the quantity of pepsin
					if( this->pepsinAlwaysDying == false )
						currentPepsin *= pepsinDyingRatio;
				}
				/*
				if( currentPosition == startingPosition && proteinsIndex == startingProtein )
				{
					attempts++;
					failedAttempts++;
					if(verbose) cout << "!!! Failed cut attempt #" << failedAttempts << " !!!" << endl;
				}
				*/
			}
		}
		else
		{
			// the pepsin did not activate
			cout << "Pepsin did not activate." << endl;
		}
	
		// if we are working with the idea that pepsin dies out with time, do it!
		if( this->pepsinAlwaysDying == true )
			currentPepsin *= pepsinDyingRatio;
		
		// increase t
		t++;
		
		// TODO momentarily removed to save memory
		// now, update the statistics by difference with the previous instant of time
		/*
		for(map< string, map<unsigned int, unsigned int> >::iterator    it = statistics.begin();
										it != statistics.end();
										it++)
		{
			// just copy the next to last quantity for each protein
			if( it->second.find( t-1 ) != it->second.end() )
			{
				unsigned int quantity = it->second[ t-1 ];
				it->second[t] = quantity;
			}
		}
		*/
		
		// if there's been a cut, modify statistics accordingly
		if( cutPerformed == true )
		{
			// first of all, increase time2
			time2++;

			// the quantity of the starting protein is lowered by 1
			lastQuantity[parentProtein.peptide] -= 1;
			statistics[ parentProtein.peptide ][t] = lastQuantity[parentProtein.peptide];

			// then, we add the two 'new' proteins
			if( lastQuantity.find( product1.peptide ) == lastQuantity.end() )
			{
				statistics[ product1.peptide ][ t ] = 1;
				lastQuantity[ product1.peptide ] = 1;
			}
			else
			{
				statistics[ product1.peptide ][ t ] = lastQuantity[product1.peptide] + 1;
				lastQuantity[product1.peptide] += 1;
			}

			if( lastQuantity.find( product2.peptide ) == lastQuantity.end() )
			{
				statistics[ product2.peptide ][ t ] = 1;
				lastQuantity[ product2.peptide ] = 1;
			}
			else
			{
				statistics[ product2.peptide ][ t ] = lastQuantity[product2.peptide] + 1;
				lastQuantity[product2.peptide] += 1;
			}
		}
		
		// take note of pepsin quantity and time2
		this->pepsinHistory.push_back( currentPepsin );
		this->time2History.push_back( time2 );
	
	} // end time loop

	return;
}

//double PepsinModel::computeCutProbability( string key, string* protein, unsigned int position )
double PepsinModel::computeCutProbability( string key, Peptide* protein, unsigned int position )
{
	if(DEBUG) cout 	<< "Now computing probability for key=\"" << key 
			<< "\", protein=" << *protein
			<< ", position=" << position 
			<< endl;
	
	// NOTE: for the index corrections in this function, you have to take into account that
	//	 "position" actually refers to P1, in classic literature; so
	//	position-1 = P2 | position+2 = P2'
	//	position-2 = P3 | position+3 = P3'
	//	position-3 = P4 | position+4 = P4'

	// base probability is 0
	double probability = 0.0;
	
	// first thing, check whether there is a disulfide bond (or similar) in the target position
	if( protein->isCuttable( position, position+1 ) == false )
	{
		if(verbose) cout 	<< "Position #" << position 
					<< " cannot be cut, due to the presence of a disulfide bond."
					<< endl;
		return probability;
	}

	// the base probability is just taken from the map
	map<string,double>::iterator it = this->cutProbability.find( key );

	// then, it is altered depending on the alterations map
	if( it != cutProbability.end() )
	{
		// compute the true probability, taking into account the alterations
		probability = it->second;

		if( alterations.size() > 0 )
		{
			if(DEBUG) cout << "- Evaluating alterations on the left..." << endl;
			// first, alterations on the left
			for(int left = 2; left <= 4 && (int)(position - left+1) >= 0; left++)
			{
				if(DEBUG) cout 	<< "- Now evaluating position " << (position - left+1) 
						<< " on protein \"" << *protein << "\""
						<< endl;

				string aminoacid = "";
				aminoacid += (*protein)[position-left+1];
				string key = "-";
				key += left + 48; // 48 is the ascii code for '0'
				
				// are there alterations for that column?
				if(DEBUG) cout 	<< "- Looking for alterations on column \"" << key 
						<< "\" for aminoacid \"" << aminoacid << "\""
						<< endl;
				map< string, map<string, double> >::iterator it2 = this->alterations.find( key );
				if( it2 != this->alterations.end() )
				{
					// are there alterations in that column, for that aminoacid?
 					map<string, double>::iterator it3 = it2->second.find( aminoacid );
					if( it3 != it2->second.end() ) probability *= alterations[key][aminoacid];
				}
			}

			// then, alterations on the right
			if(DEBUG) cout << "- Evaluating alterations on the right..." << endl;
			for(unsigned int right = 2; right <= 4 && 
				(position + right) < protein->length(); right++)
			{
				string aminoacid = "";
				aminoacid += (*protein)[position+right];
				string key = "+";
				key += right + 48; // 48 is the ascii code for '0'

				// are there alterations for that column?
				map< string, map<string, double> >::iterator it2 = this->alterations.find( key );
				if( it2 != this->alterations.end() )
				{
					// are there alterations in that column, for that aminoacid?
 					map<string, double>::iterator it3 = it2->second.find( aminoacid );
					if( it3 != it2->second.end() ) probability *= alterations[key][aminoacid];
				}
			}
			
			// finally, let's check whether there is a missing peptide (meaning that we are either
			// at the beginning, or at the end of a chain
			if(DEBUG) cout << "- Evaluating alterations close to an end of the chain..." << endl;
			double multiplier = 1.0;
			
			// iterate over the keys
			for(	std::map<int, double>::iterator terminalIt = this->terminalAlterations.begin();
				terminalIt != this->terminalAlterations.end(); terminalIt++
			)
			{
				// check if we are at the end of a chain
				int key = terminalIt->first;
				double keyMultiplier = terminalIt->second;
				// left part (key < 0)	      right part (key > 0)
				if( position + key + 1 < 0 || position + key >= protein->length() )
				{
					if(DEBUG) cout 	<< "For position=" << position << " (" << protein->peptide[position]
							<< ") and key=" << key
							<< " in protein=\"" << protein->peptide << "\", we found an end!"
							<< endl;
					
					// only the LOWEST multiplier applies
					if( keyMultiplier < multiplier ) multiplier = keyMultiplier;
				}
			}
			
			// apply the multiplier
			probability *= multiplier;
			

		} // end if alterations size > 0
	} // end if there are corresponding probabilities in the table
	
	return probability;
}
