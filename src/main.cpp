// Main for the model

// standard library
#include <iostream>
#include <cstdio>
#include <cstring>
#include <vector>
#include <memory>

// custom classes
#include <PepsinModel.h>
#include <Log.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

// defines
#define OP_INPUT "--input"
#define OP_OUTPUT "--output"
#define OP_REPETITIONS "--repetitions"
#define OP_VERBOSE "--verbose"
#define OP_LOGFILE "--log-file"
#define OP_LOGLEVEL "--log-level"

#define DEFAULT_OUTPUT "statistics.csv"
#define DEFAULT_REPETITIONS 1

#define EXIT_SUCCESS 0

using namespace std;

void printUsage(char* programName)
{
	cout 	<< endl << "Usage: " << programName << " "
		<< OP_INPUT << " <input.json> "
		<< "[" << OP_OUTPUT << " <statistics.csv> "
		<< OP_VERBOSE << " "
		<< OP_LOGFILE << " <debug.log> "
		<< OP_LOGLEVEL << " <level> ]" << endl
		<< endl
		<< "Input arguments:" << endl
		<< "\t" << OP_INPUT << " <input.json>: the JSON input file with the complete description of the problem." << endl
		<< "\t" << OP_OUTPUT << " <statistics.csv>: the CSV output file with the solution of the problem. This argument is optional, the default name is \"" << DEFAULT_OUTPUT << "\"" << endl
		<< "\t" << OP_VERBOSE << " : optional argument, equivalent to \"" << OP_LOGLEVEL << " debug\" on the console. By default, only info/warning/error messages are shown." << endl
		<< "\t" << OP_LOGFILE << " <debug.log> : optional argument, also writes a full trace-level debug log to the given file. No file is created unless this is specified." << endl
		<< "\t" << OP_LOGLEVEL << " <level> : optional argument, sets the console log level explicitly (one of: trace, debug, info, warn, error, off). Overrides " << OP_VERBOSE << " if both are given." << endl
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
	string logFile;
	string logLevelArg;

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
		else if( strcmp(argv[a], OP_LOGFILE) == 0 && a+1 < argc )
		{
			logFile = argv[a+1];
		}
		else if( strcmp(argv[a], OP_LOGLEVEL) == 0 && a+1 < argc )
		{
			logLevelArg = argv[a+1];
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

	// figure out the console level: --log-level wins if given, otherwise --verbose means
	// "debug", and the default is "info" (quiet: only progress/warnings/errors, no trace spam)
	spdlog::level::level_enum consoleLevel = verbose ? spdlog::level::debug : spdlog::level::info;
	if( !logLevelArg.empty() )
	{
		spdlog::level::level_enum parsedLevel = spdlog::level::from_str( logLevelArg );
		if( parsedLevel == spdlog::level::off && logLevelArg != "off" )
		{
			cerr 	<< "Warning: log level \"" << logLevelArg << "\" not recognized "
				<< "(valid values are: trace, debug, info, warn, error, off). Using the default instead." << endl;
		}
		else
		{
			consoleLevel = parsedLevel;
		}
	}

	// console sink is always present; a file sink is only added if --log-file was specified,
	// so no text file is ever written to the user's folder unless explicitly asked for
	vector<spdlog::sink_ptr> sinks;

	auto consoleSink = make_shared<spdlog::sinks::stdout_color_sink_mt>();
	consoleSink->set_level( consoleLevel );
	sinks.push_back( consoleSink );

	// the logger's own level must be at least as permissive as the most verbose sink;
	// each sink then filters independently via its own set_level() above
	spdlog::level::level_enum loggerLevel = consoleLevel;

	if( !logFile.empty() )
	{
		auto fileSink = make_shared<spdlog::sinks::basic_file_sink_mt>( logFile, /*truncate*/ true );
		fileSink->set_level( spdlog::level::trace );
		sinks.push_back( fileSink );
		loggerLevel = spdlog::level::trace;
	}

	logging::setSinks( sinks );
	logging::setLevel( loggerLevel );

	// create an instance of the model
	PepsinModel pepsinModel;

	// load JSON configuration
	if( pepsinModel.readJson( jsonInputFile ) != EXIT_SUCCESS )
	{
		return -1;
	}

	// run the model
	pepsinModel.run();

	// write statistics to file
	LOG_INFO("Writing statistics to file \"" << csvOutputFile << "\"...");
	pepsinModel.writeLog( csvOutputFile );

	// end
	LOG_INFO("Done.");

	return 0;
}
