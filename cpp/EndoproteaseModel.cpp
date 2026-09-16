// Methods for the class representing a generic endoprotease model

// class header
#include <EndoproteaseModel.h>

// other classes used by this model
#include <Peptide.h>

// standard classes
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <ctime>
#include <sstream>
#include <vector>
#include <random>

// local classes/libraries
// (JSON parsing; used to be tinyxml, until the switch from XML to JSON configuration files)
#include <nlohmann/json.hpp>
using json = nlohmann::json;

// leveled logging (replaces the old "verbose" bool and hardcoded "DEBUG" macro)
#include <Log.h>

// some defines, used mainly for the JSON parsing
#define DEFAULT_CSV "statistics-"
#define DEFAULT_OUTPUT "solution.xml"
#define DEFAULT_REPETITIONS 1

#define MAX_ATTEMPTS 10
#define MAX_FAILED_ATTEMPTS 100
#define MAX_TIME 10000

#define NAME_ENZYME "enzyme"

#define JSON_ALTERATIONS "alterations"
#define JSON_CUTS "cuts"
#define JSON_INITIALENZYME "initialEnzyme"
#define JSON_MAXATTEMPTS "maxAttempts"
#define JSON_MAXATTEMPTSPERTIME "maxAttemptsPerTime"
#define JSON_MAXDH "maxDH"
#define JSON_MAXTIME "maxTime"
#define JSON_PARAMETERS "parameters"
#define JSON_ENZYMEALWAYSDYING "enzymeAlwaysDying"
#define JSON_ENZYMEDYINGRATIO "enzymeDyingRatio"
#define JSON_PROTEINS "proteins"
#define JSON_RANDOMSEED "randomSeed"
#define JSON_TERMINALALTERATIONS "terminalAlterations"

#define JSON_FIELD_DISULFIDEBONDS "disulfideBonds"
#define JSON_FIELD_LEFT "left"
#define JSON_FIELD_NAME "name"
#define JSON_FIELD_QUANTITY "quantity"
#define JSON_FIELD_RIGHT "right"
#define JSON_FIELD_SEQUENCE "sequence"

// namespace
using namespace std;

// constructor (basically it does nothing at all, besides providing a few default values)
EndoproteaseModel::EndoproteaseModel() :
currentEnzyme(1.0),
maxAttemptsPerTime(10),
maxAttempts(1000),
maxDH(100.0),
maxTime(10000),
overallLength(0),
enzymeAlwaysDying(false),
enzymeDyingRatio(1.0),
randomSeed(0),
t(0)
{}

// destructor
EndoproteaseModel::~EndoproteaseModel()
{
	// clear all the maps and vectors
	this->alterations.clear();
	this->cutProbability.clear();
	this->originalProteins.clear();
	this->enzymeHistory.clear();
	this->statistics.clear();
}

// draw a random double in [0, 1) from this instance's own random engine
//
// NOTE: deliberately not std::uniform_real_distribution. std::mt19937 itself is fully
// specified by the standard (same output, same seed, on every platform/compiler), but the
// *distributions* built on top of an engine are only implementation-defined -- libstdc++
// (Linux, MinGW) and libc++ (macOS/Clang) map the same engine output to different final
// values. Found via CI: a fixed-seed regression test matched on Linux (MinGW and Linux GCC
// both use libstdc++) but produced a different-shaped result on macOS. Computing the [0,1)
// mapping by hand, directly from the engine's own (portable) min()/max()/operator()(), keeps
// this fully portable across standard library implementations too, not just across
// platforms that happen to share one.
double EndoproteaseModel::randomUnit()
{
	double range = static_cast<double>( randomEngine.max() - randomEngine.min() ) + 1.0;
	return static_cast<double>( randomEngine() - randomEngine.min() ) / range;
}

// draw a random unsigned int in [0, exclusiveUpperBound) from this instance's own random
// engine; see randomUnit() above for why this doesn't use std::uniform_int_distribution
unsigned int EndoproteaseModel::randomIndex( unsigned int exclusiveUpperBound )
{
	return static_cast<unsigned int>( randomUnit() * exclusiveUpperBound );
}

// read a JSON configuration file from disk, then hand off to readJsonObject() below.
// Kept separate so the same parsing logic can also be driven directly from an
// already-in-memory JSON value (e.g. a Python dict, converted via pybind11_json) without
// needing to write it to a temporary file first.
int EndoproteaseModel::readJson( string fileName )
{
	LOG_DEBUG("Loading JSON file...");

	ifstream inStream( fileName.c_str() );
	if( !inStream.is_open() )
	{
		LOG_ERROR("Error: cannot open file \"" << fileName << "\". Aborting...");
		return -1;
	}

	json root;
	try
	{
		// ignore_comments=true so the config file can use "//" and "/* */" comments,
		// just like the old XML file was heavily commented
		root = json::parse( inStream, /*callback*/ nullptr, /*allow_exceptions*/ true, /*ignore_comments*/ true );
	}
	catch( json::parse_error& e )
	{
		LOG_ERROR("Error: file \"" << fileName << "\" is not valid JSON (" << e.what() << "). Aborting...");
		return -1;
	}

	return this->readJsonObject( root );
}

// store the relevant information from an already-parsed JSON configuration
int EndoproteaseModel::readJsonObject( const json& root )
{
	// first of all, parse the "parameters"
	// parameters are not compulsory (TODO REALLY?)
	if( root.contains(JSON_PARAMETERS) )
	{
		const json& parameters = root[JSON_PARAMETERS];

		// check if there's a random seed (null / absent means "use time-based seed")
		if( parameters.contains(JSON_RANDOMSEED) && !parameters[JSON_RANDOMSEED].is_null() )
			this->randomSeed = parameters[JSON_RANDOMSEED].get<unsigned int>();

		if( parameters.contains(JSON_MAXTIME) )
			this->maxTime = parameters[JSON_MAXTIME].get<unsigned int>();

		if( parameters.contains(JSON_MAXDH) )
			this->maxDH = parameters[JSON_MAXDH].get<double>();

		if( parameters.contains(JSON_MAXATTEMPTSPERTIME) )
			this->maxAttemptsPerTime = parameters[JSON_MAXATTEMPTSPERTIME].get<unsigned int>();

		if( parameters.contains(JSON_MAXATTEMPTS) )
			this->maxAttempts = parameters[JSON_MAXATTEMPTS].get<unsigned int>();

		if( parameters.contains(JSON_INITIALENZYME) )
			this->currentEnzyme = parameters[JSON_INITIALENZYME].get<double>();

		// check if the enzyme dies at every timestep, or just when it does not cut
		// (NOTE: the old XML parser had a copy-paste bug here, reading this value from the
		// wrong field entirely; fixed as part of the JSON migration)
		if( parameters.contains(JSON_ENZYMEALWAYSDYING) )
			this->enzymeAlwaysDying = parameters[JSON_ENZYMEALWAYSDYING].get<bool>();

		if( parameters.contains(JSON_ENZYMEDYINGRATIO) )
			this->enzymeDyingRatio = parameters[JSON_ENZYMEDYINGRATIO].get<double>();
	}

	if( randomSeed != 0 )
	{
		// initialize this instance's random engine with the given seed
		LOG_DEBUG("Initializing random generator with seed " << randomSeed);
		this->randomEngine.seed( randomSeed );
	}
	else
	{
		// initialize this instance's random engine with the current time
		time_t timeRandomSeed = time(NULL);
		LOG_DEBUG("Initializing random generator with time (" << timeRandomSeed << ")");
		this->randomEngine.seed( static_cast<unsigned int>(timeRandomSeed) );
	}

	// take the proteins
	// include the possibility of having different proteins in the initial set
	if( !root.contains(JSON_PROTEINS) )
	{
		LOG_ERROR("Error: field \"" << JSON_PROTEINS << "\" must be specified. Aborting...");
		return -1;
	}

	for( const json& protein : root[JSON_PROTEINS] )
	{
		LOG_DEBUG("Adding a protein...");
		string originalProtein = protein.value(JSON_FIELD_SEQUENCE, string(""));

		// ok, here we get two possible fields: the name of the protein (unused) and the quantity (important!)
		// recently, we added a third field, a list of disulfide bonds
		if( protein.contains(JSON_FIELD_NAME) )
			LOG_DEBUG("Adding protein \"" << protein[JSON_FIELD_NAME].get<string>() << "\"...");

		// default quantity for each protein is 1
		unsigned int proteinQuantity = protein.value(JSON_FIELD_QUANTITY, 1u);

		// remove all whitespaces and other stuff from the string
		originalProtein.erase( std::remove_if( originalProtein.begin(), originalProtein.end(), ::isspace ), originalProtein.end() );

		LOG_DEBUG(	"Adding protein with composition \"" << originalProtein
				<< "\", " << originalProtein.length()
				<< " amminoacids long." );

		// create the Peptide with the original protein
		Peptide originalProteinPeptide(originalProtein);

		// add disulfide bonds (if they're there); positions in the file are 1-indexed
		if( protein.contains(JSON_FIELD_DISULFIDEBONDS) )
		{
			for( unsigned int bond : protein[JSON_FIELD_DISULFIDEBONDS] )
			{
				LOG_TRACE(	"Adding bond in position #" << bond
						<< " to the protein (modified to #"
						<< (bond-1) << " for C++ internal array indexing)" );

				originalProteinPeptide.addDisulfideBond( bond-1 );
			}
		}

		for(unsigned int i = 0; i < proteinQuantity; i++)
		{
			LOG_TRACE("\tAdding copy #" << (i+1) << "...");
			// add the protein to the initial set
			this->originalProteins.push_back( originalProteinPeptide );

			// increase the overall length of the original proteins
			this->overallLength += originalProtein.length();
		}
	}

	// check if there are proteins at all
	if( originalProteins.size() == 0 )
	{
		LOG_ERROR("Error: no proteins found in \"" << JSON_PROTEINS << "\". Aborting...");
		return -1;
	}

	// then, take all the cuts: "cuts": { "<left>": { "<right>": probability, ... }, ... }
	LOG_DEBUG("Reading the cuts...");
	if( !root.contains(JSON_CUTS) )
	{
		LOG_ERROR("Error: field \"" << JSON_CUTS << "\" must be specified. Aborting...");
		return -1;
	}

	for( auto& leftEntry : root[JSON_CUTS].items() )
	{
		const string& left = leftEntry.key();
		for( auto& rightEntry : leftEntry.value().items() )
		{
			const string& right = rightEntry.key();
			double probability = rightEntry.value().get<double>();

			// put everything inside the map, if probability > 0
			if( probability > 0.0 )
			{
				string key = left;
				key += "\t";
				key += right;
				this->cutProbability[ key ] = probability;

				LOG_TRACE("Cut \"" << key << "\" with probability=" << probability << " added");
			}
		}
	}

	// check if the map is empty, it could be a problem
	if( cutProbability.size() == 0 )
	{
		LOG_WARN("Warning: there are no cuts specified, so probably nothing will happen during this simulation...");
	}
	else
	{
		LOG_DEBUG("There are " << cutProbability.size() << " cuts specified");
	}

	// proceed with managing the alterations
	// alterations are optional, so if they're not there it's not a big deal
	if( root.contains(JSON_ALTERATIONS) )
	{
		LOG_DEBUG("Now reading alterations...");

		// "alterations": { "<aminoacid>": { "left": { "<position>": probability, ... }, "right": {...} }, ... }
		for( auto& aminoacidEntry : root[JSON_ALTERATIONS].items() )
		{
			const string& aminoacid = aminoacidEntry.key();
			const json& sides = aminoacidEntry.value();

			for( const string& side : { string(JSON_FIELD_LEFT), string(JSON_FIELD_RIGHT) } )
			{
				if( !sides.contains(side) ) continue;

				for( auto& positionEntry : sides[side].items() )
				{
					string key = (side == JSON_FIELD_LEFT) ? "-" : "+";
					key += positionEntry.key();
					double probability = positionEntry.value().get<double>();

					// store the alteration inside the map
					this->alterations[ key ][ aminoacid ] = probability;
					LOG_TRACE(	"- alteration[ " << key << " ][ " << aminoacid << " ] = "
							<< this->alterations[key][aminoacid] );
				}
			}
		}

		// we also have a second type of alteration, called "terminal", to take into account
		// the alteration of the probabilities of cutting, if you are at either end of a chain
		// "terminalAlterations": { "left": { "<position>": multiplier, ... }, "right": {...} }
		LOG_DEBUG("Now reading terminal alterations...");
		if( root.contains(JSON_TERMINALALTERATIONS) )
		{
			const json& terminalAlterationsJson = root[JSON_TERMINALALTERATIONS];

			for( const string& side : { string(JSON_FIELD_LEFT), string(JSON_FIELD_RIGHT) } )
			{
				if( !terminalAlterationsJson.contains(side) ) continue;

				for( auto& positionEntry : terminalAlterationsJson[side].items() )
				{
					int key = std::stoi( positionEntry.key() );
					double multiplier = positionEntry.value().get<double>();

					if( side == JSON_FIELD_LEFT )
						key = -1 * key;

					// store this alteration inside the map
					this->terminalAlterations[ key ] = multiplier;
					LOG_DEBUG("- terminal alteration[ " << key << " ] = " << this->terminalAlterations[key]);
				}
			}
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

			LOG_DEBUG(	"For alterations on column \"" << it->first
					<< "\", avg=" << average
					<< ", std=" << std
					<< ", sigma=" << sqrt(std) );

			// remake the hash map, inserting only the values that will actually modify the probabilities;
			// the alterations might be
			for( map<string, double>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++ )
			{
				double value = it2->second - average;
				if( value < 0 ) value = -1.0 * value;

				if( value > std )
				{
					LOG_DEBUG(	"- alteration[ " << it->first << " ][ "
							<< it2->first << " ]=" << it2->second
							<< ", difference from average is " << (value / std ) << " variances!" );
				}

				if( value > sqrt(std) )
				{
					LOG_DEBUG(	"- alteration[ " << it->first << " ][ "
							<< it2->first << " ]=" << it2->second
							<< ", difference from average is " << (value / sqrt(std) ) << " standard deviations!" );
				}

				it2->second /= average;
			}
		}


		// this is just debugging / pre-processing, and it should probably be removed; but it is interesting
		// to understand how the cut probabilities are changed by the alterations: so we iterate over the
		// proteins and check the difference with the proposed approach
		if( logging::get()->should_log(spdlog::level::trace) )
		for(unsigned int p = 0; p < 1 /* originalProteins.size() */; p++)
		{
			LOG_TRACE("Analyzing the influence of alterations for protein \"" << originalProteins[p] << "\"...");
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
					LOG_TRACE("- \"" << key << "\": original cut probability = " << probability);

					probability = this->computeCutProbability( key, &originalProteins[p], i );
					LOG_TRACE("- \"" << key << "\": altered cut probability = " << probability << "\n");
				}

			}
		}
		
		
	} // end if there are alterations
	
	
	// everything all right!
	return 0;
}

// Build the (row, column) time series used for output: one row per "print point" (every
// `period`-th distinct value of time2, to avoid HUGE outputs -- see the note below), one
// column per distinct peptide produced during the simulation, counts carried forward between
// change-points.
//
// This used to be folded directly into writeLog(), re-checking (per-column, per-localt)
// whether the current iteration was a print point -- for the full lactoferrin run that meant
// something like 990,000 iterations x 10,000 peptide columns = ~9.9 BILLION map::find() calls
// to produce just 3,441 output rows, which is what made that function take ~15 minutes.
// TODO
// - is it possible to change the increments so that localt += period?
// - period readable by configuration file
EndoproteaseModel::TimeSeries EndoproteaseModel::computeTimeSeries( unsigned int period )
{
	LOG_DEBUG("Pre-processing statistics...");

	TimeSeries series;

	// figure out, in a single cheap pass over 0..t, exactly which localt values are actually
	// going to be included: the first localt to reach each new multiple of `period` in time2
	{
		int lastTime2 = -1;
		for(unsigned int localt = 0; localt < t; localt++)
		{
			if( this->time2History[localt] % period == 0 && (int)this->time2History[localt] != lastTime2 )
			{
				series.time.push_back( localt );
				series.time2.push_back( this->time2History[localt] );
				series.enzyme.push_back( this->enzymeHistory[localt] );
				lastTime2 = this->time2History[localt];
			}
		}
	}

	// for each peptide column, walk its own (sparse) map of change-points alongside series.time
	// with a simple iterator advance instead of a fresh find() per row: both are sorted by t,
	// so this is a linear merge, and each column ends up doing roughly (series.time.size() +
	// its own number of changes) work in total, instead of t.size() lookups regardless of
	// whether the row is even going to be included
	series.peptideNames.reserve( statistics.size() );
	series.peptideCounts.assign( statistics.size(), vector<unsigned int>( series.time.size(), 0 ) );

	size_t col = 0;
	for(map<string, map<unsigned int, unsigned int> >::iterator 	it = statistics.begin();
									it != statistics.end();
									it++, col++)
	{
		series.peptideNames.push_back( it->first );

		unsigned int columnValue = 0;
		map<unsigned int, unsigned int>::const_iterator changeIt = it->second.begin();
		map<unsigned int, unsigned int>::const_iterator changeEnd = it->second.end();

		for(size_t row = 0; row < series.time.size(); row++)
		{
			unsigned int localt = series.time[row];

			// advance to the most recent change at or before this row's t
			while( changeIt != changeEnd && changeIt->first <= localt )
			{
				columnValue = changeIt->second;
				++changeIt;
			}

			series.peptideCounts[col][row] = columnValue;
		}
	}

	return series;
}

// write the history to a CSV file
int EndoproteaseModel::writeLog( string fileName )
{
	TimeSeries series = this->computeTimeSeries();

	// trying to put some buffering, in order not to block the file for too long
	stringstream outStream;

	// header
	outStream << "\"time\",\"time2\",\"" << NAME_ENZYME << "\"";
	for(const string& name : series.peptideNames)
		outStream << ",\"" << name << "\"";
	outStream << "\n";

	// data rows
	for(size_t row = 0; row < series.time.size(); row++)
	{
		outStream 	<< series.time[row] << ","
				<< series.time2[row] << ","
				<< series.enzyme[row];

		for(size_t col = 0; col < series.peptideNames.size(); col++)
			outStream << "," << series.peptideCounts[col][row];

		outStream << "\n";
	}

	LOG_DEBUG("Writing statistics to CSV file \"" << fileName << "\"...");
	// open file
	ofstream csvOut( fileName.c_str() );
	if( !csvOut.is_open() )
	{
		LOG_ERROR("Error: cannot write on file \"" << fileName << "\". Aborting...");
		return -1;
	}

	csvOut << outStream.str();

	csvOut.close();

	// everything all right!
	return 0;
}

// after initialization, here is the "true" run of the model!
void EndoproteaseModel::run()
{
	LOG_DEBUG(	"Starting the simulation, with maxTime=" << maxTime
			<< ", maxAttemptsPerTime=" << maxAttemptsPerTime
			<< ", maxAttempts=" << maxAttempts
			<< ", maxDH=" << maxDH
			<< ", initialEnzyme=" << currentEnzyme
			<< ", enzymeDyingRatio=" << enzymeDyingRatio );

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
		// the total number of sites where the enzyme can cut the proteins
		dhTotalSites += (uniqueProteins[i].length() - 1) * occurrences;
	}

	// some debugging
	LOG_DEBUG("The total number of sites that the enzyme can cut is " << dhTotalSites);

	// also, take note of enzyme quantity and time2
	this->enzymeHistory.push_back( currentEnzyme );
	this->time2History.push_back( time2 );

	// start!
	while( this->t < this->maxTime && attempts < this->maxAttempts && (time2/(double)dhTotalSites) < this->maxDH )
	{
		LOG_DEBUG(	"Time #" << t
				<< " (Time2 #" << time2
				<< ", DH=" << (time2/(double)dhTotalSites)
				<< ", maxDH=" << this->maxDH
				<< ")" );
		
		// some variables need to be defined here
		bool cutPerformed = false;
		//string parentProtein, product1, product2;
		Peptide parentProtein, product1, product2;

		// reset the attempts done per instant of t
		unsigned int attemptsPerTime = 0;

		// before doing all the rest, just verify whether the enzyme acts
		double randomEnzymeActivation = this->randomUnit();

		if( randomEnzymeActivation < currentEnzyme)
		{
			// find a random point in a random protein
			// TODO this part could be sped up by creating a map< unsigned int, pair<unsigned int, unsigned int> > 
			//	that keeps track of the protein and relative position for an absolute position;
			//	so, if I have 200 peptides with 20 aminoacids each, I have a total of 4000 positions;
			//	the map can tell me that position 201 is actually position 0 on protein 1

			// ok, *probably* the probability of a random protein to be chosen is 
			// dependent on its size; so, we should probably just choose a random point in the starting protein and
			// go over all the proteins until you find the corresponding point
			unsigned int startingPosition = this->randomIndex( overallLength );
			
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
			
			LOG_DEBUG(	"Starting position is " << startingPosition
					<< " in protein #" << startingProtein
					<< " (\"" << currentProteins[startingProtein] << "\")" );

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
				
				LOG_TRACE(	"- Analyzing protein[" << proteinsIndex << "][" << currentPosition
						<< "]=\"" << key << "\"..." );

				// obtain the probability of cutting in that point
				double probability = this->computeCutProbability( key, &currentProteins[proteinsIndex], currentPosition );

				if( probability > 0.0 )
				{
					// random number between 0 and 1
					double randomNumber = this->randomUnit();
						
					if( randomNumber < probability )
					{
						// perform the cut!
						LOG_DEBUG(	"Performing cut \"" << key << "\" in position "
								<< currentPosition << " that had probability "
								<< cutProbability[key]
								<< " (modified to " << probability
								<< ", randomNumber was " << randomNumber << ")" );
						
						// now, the current proteins will change! add the resulting sub-proteins
						// at the end of the vector, and remove the current protein
						product1 = currentProteins[proteinsIndex].substr(0, currentPosition+1);
						product2 = currentProteins[proteinsIndex].substr(currentPosition+1);
						
						LOG_DEBUG(	"Protein \"" << currentProteins[proteinsIndex]
								<< "\" has been cut in position " << currentPosition
								<< " and the two resulting products are \"" << product1
								<< "\" and \"" << product2 << "\"" );

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
						LOG_TRACE(	"Probability was " << probability
								<< ", while randomNumber was " << randomNumber
								<< ": cut not performed." );
					}

				}
				else
				{
					LOG_TRACE("Probability to cut = 0.0; there are no cuts for this key, or there's a disulfide bond near the target link.");
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
						LOG_DEBUG(	"Going over the list of proteins again: attemptsPerTime="
								<< attemptsPerTime << ", attempts=" << attempts
								<< ", startingPosition=" << startingPosition
								<< ", startingProtein=" << startingProtein );
					}


					currentPosition = 0;
					LOG_TRACE(	"- Now analyzing protein #" << proteinsIndex
							<< " (size " << currentProteins[proteinsIndex].length() << ")" );
				}

				// if we came back where we started...
				if( cutPerformed == false )
				{
					attemptsPerTime++;
					attempts++;
					
					// also, reduce the quantity of enzyme
					if( this->enzymeAlwaysDying == false )
						currentEnzyme *= enzymeDyingRatio;
				}
				/*
				if( currentPosition == startingPosition && proteinsIndex == startingProtein )
				{
					attempts++;
					failedAttempts++;
					LOG_TRACE("!!! Failed cut attempt #" << failedAttempts << " !!!");
				}
				*/
			}
		}
		else
		{
			// the enzyme did not activate
			LOG_DEBUG("Enzyme did not activate.");
		}

		// if we are working with the idea that the enzyme dies out with time, do it!
		if( this->enzymeAlwaysDying == true )
			currentEnzyme *= enzymeDyingRatio;
		
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
		
		// take note of enzyme quantity and time2
		this->enzymeHistory.push_back( currentEnzyme );
		this->time2History.push_back( time2 );
	
	} // end time loop

	return;
}

//double EndoproteaseModel::computeCutProbability( string key, string* protein, unsigned int position )
double EndoproteaseModel::computeCutProbability( string key, Peptide* protein, unsigned int position )
{
	LOG_TRACE(	"Now computing probability for key=\"" << key
			<< "\", protein=" << *protein
			<< ", position=" << position );
	
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
		LOG_TRACE(	"Position #" << position
				<< " cannot be cut, due to the presence of a disulfide bond." );
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
			LOG_TRACE("- Evaluating alterations on the left...");
			// first, alterations on the left
			for(int left = 2; left <= 4 && (int)(position - left+1) >= 0; left++)
			{
				LOG_TRACE(	"- Now evaluating position " << (position - left+1)
						<< " on protein \"" << *protein << "\"" );

				string aminoacid = "";
				aminoacid += (*protein)[position-left+1];
				string key = "-";
				key += left + 48; // 48 is the ascii code for '0'

				// are there alterations for that column?
				LOG_TRACE(	"- Looking for alterations on column \"" << key
						<< "\" for aminoacid \"" << aminoacid << "\"" );
				map< string, map<string, double> >::iterator it2 = this->alterations.find( key );
				if( it2 != this->alterations.end() )
				{
					// are there alterations in that column, for that aminoacid?
 					map<string, double>::iterator it3 = it2->second.find( aminoacid );
					if( it3 != it2->second.end() ) probability *= alterations[key][aminoacid];
				}
			}

			// then, alterations on the right
			LOG_TRACE("- Evaluating alterations on the right...");
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
			LOG_TRACE("- Evaluating alterations close to an end of the chain...");
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
					LOG_TRACE(	"For position=" << position << " (" << protein->peptide[position]
							<< ") and key=" << key
							<< " in protein=\"" << protein->peptide << "\", we found an end!" );
					
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
