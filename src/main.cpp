// Main for the model

// standard library
#include <iostream>
#include <cstdio>
#include <cstring>

// custom classes
#include <PepsinModel.h>

// defines
#define OP_INPUT "--input"
#define OP_OUTPUT "--output"
#define OP_REPETITIONS "--repetitions"
#define OP_VERBOSE "--verbose"

#define DEFAULT_OUTPUT "statistics.csv"
#define DEFAULT_REPETITIONS 1

#define EXIT_SUCCESS 0

using namespace std;

void printUsage(char* programName)
{
	cout 	<< endl << "Usage: " << programName << " "
		<< OP_INPUT << " <input.json> "
		<< "[" << OP_OUTPUT << " <statistics.csv> "
		<< OP_VERBOSE << " ]" << endl
		<< endl
		<< "Input arguments:" << endl
		<< "\t" << OP_INPUT << " <input.json>: the JSON input file with the complete description of the problem." << endl
		<< "\t" << OP_OUTPUT << " <statistics.csv>: the CSV output file with the solution of the problem. This argument is optional, the default name is \"" << DEFAULT_OUTPUT << "\"" << endl
		<< "\t" << OP_VERBOSE << " : optional argument, regulates the verbosity of the program. By default, it's silent if not for possible errors." << endl
		<< endl;

	return;
}

// main function
int main(int argc, char* argv[])
{
	// some variables
	bool verbose = false;
	unsigned int repetitions = DEFAULT_REPETITIONS;

	string jsonInputFile;
	string csvOutputFile;

	// parse arguments
	for(unsigned int a = 0; a < argc; a++)
	{
		if( strcmp(argv[a], OP_INPUT) == 0 && a+1 < argc )
		{
			jsonInputFile = argv[a+1];
		}
		else if( strcmp(argv[a], OP_OUTPUT) == 0 && a+1 < argc )
		{
			csvOutputFile = argv[a+1];
		}
		else if( strcmp(argv[a], OP_REPETITIONS) == 0 && a+1 < argc )
		{
			if( sscanf( argv[a+1], "%u", &repetitions ) != 1 )
			{
				cerr 	<< "Warning: value \"" << argv[a+1] << "\" specified for the repetitions is not valid. "
					<< "The default value (" << repetitions << ") will be used instead." << endl;
			}
		}
		else if( strcmp(argv[a], OP_VERBOSE) == 0 )
		{
			verbose = true;
		}
		else if( argv[a][0] == '-' )
		{
			cerr << "Warning: option \"" << argv[a] << "\" not recognized..." << endl;
		}
	}

	// inputs check
	if( jsonInputFile.length() <= 0 )
	{
		cerr << "Error: a JSON input file must be specified. Aborting..." << endl;
		printUsage( argv[0] );
		return 0;
	}
	
	if( csvOutputFile.length() <= 0 )
	{
		csvOutputFile = DEFAULT_OUTPUT;
	}
	
	// create an instance of the model
	PepsinModel pepsinModel;
	
	// set verbosity of the model
	pepsinModel.verbose = verbose;
	
	// load JSON configuration
	if( pepsinModel.readJson( jsonInputFile ) != EXIT_SUCCESS )
	{
		return -1;
	}
	
	// run the model
	pepsinModel.run();

	// write statistics to file
	cout << "Writing statistics to file \"" << csvOutputFile << "\"..." << endl;
	pepsinModel.writeLog( csvOutputFile );

	// end
	cout << "Done." << endl;
	
	return 0;
}
