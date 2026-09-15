// pybind11 bindings for the endocleave Python package.
//
// Two design points worth explaining:
//
// 1. Logging is bridged into Python's own `logging` module (see PythonLoggingSink below)
//    instead of a console/file sink, so `import endocleave` doesn't print anything on its
//    own, and Python users configure it the normal Python way. BUT: the C++ side's log
//    level is still the real performance gate (see Log.h's should_log() check, used by
//    every LOG_TRACE(...)/LOG_DEBUG(...) call site) -- it has to stay in sync with whatever
//    level is actually wanted, rather than being left permissive, otherwise every one of the
//    millions of LOG_TRACE calls in a real simulation would build its message string and
//    cross into Python (acquiring the GIL each time) only to be dropped by Python's own
//    level check. That's exactly the kind of "do the work anyway, only to throw it away"
//    bug this project already paid for once (see the writeLog() history). set_log_level()
//    below keeps both sides in sync in one call.
//
// 2. EndoproteaseModel::run() releases the GIL for its (potentially long) duration, so it
//    doesn't block other Python threads; the logging sink re-acquires the GIL just for the
//    brief moment it needs to call into Python for each message that actually passes the
//    level check above.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11_json/pybind11_json.hpp>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/base_sink.h>

#include <EndoproteaseModel.h>
#include <Log.h>

#include <mutex>
#include <string>

namespace py = pybind11;

namespace
{
	// spdlog level <-> Python's logging module level numbers (avoids importing the
	// "logging" module just to read its int constants, which are stable and documented)
	int spdlogLevelToPython( spdlog::level::level_enum level )
	{
		switch( level )
		{
			case spdlog::level::trace:
			case spdlog::level::debug:    return 10; // logging.DEBUG
			case spdlog::level::info:     return 20; // logging.INFO
			case spdlog::level::warn:     return 30; // logging.WARNING
			case spdlog::level::err:      return 40; // logging.ERROR
			case spdlog::level::critical: return 50; // logging.CRITICAL
			default:                      return 0;  // logging.NOTSET
		}
	}

	// forwards every spdlog record that passes the level check (see should_log() in Log.h)
	// into Python's logging.getLogger("endocleave")
	class PythonLoggingSink : public spdlog::sinks::base_sink<std::mutex>
	{
	public:
		PythonLoggingSink()
		{
			py::gil_scoped_acquire gil;
			pyLogger = py::module_::import("logging").attr("getLogger")( logging::LOGGER_NAME );
		}

	protected:
		void sink_it_( const spdlog::details::log_msg& msg ) override
		{
			py::gil_scoped_acquire gil;

			int pyLevel = spdlogLevelToPython( msg.level );
			std::string payload( msg.payload.data(), msg.payload.size() );

			// Logger.log(level, msg) -- no extra "args", so no %-style formatting is
			// attempted on payload even if it happens to contain a literal '%'
			pyLogger.attr("log")( pyLevel, payload );
		}

		void flush_() override {}

	private:
		py::object pyLogger;
	};

	// keeps the C++-side performance gate and Python's own logger level in sync; see the
	// file-level comment above for why both need to move together
	void setLogLevel( const std::string& levelName )
	{
		spdlog::level::level_enum level = spdlog::level::from_str( levelName );

		logging::setLevel( level );

		py::object loggingModule = py::module_::import("logging");
		loggingModule.attr("getLogger")( logging::LOGGER_NAME ).attr("setLevel")( spdlogLevelToPython(level) );
	}
}

PYBIND11_MODULE(_core, m)
{
	m.doc() = "endocleave native extension: pybind11 bindings for the C++ endoprotease simulation core";

	auto pySink = std::make_shared<PythonLoggingSink>();
	logging::setSinks( { pySink } );

	// matches the CLI's own default verbosity; call endocleave.set_log_level(...) to change it
	setLogLevel( "info" );

	m.def( "set_log_level", &setLogLevel, py::arg("level"),
		"Set the log level (\"trace\", \"debug\", \"info\", \"warning\", \"error\", or \"off\"), "
		"for both the native core and Python's logging.getLogger(\"endocleave\"). Use this "
		"instead of calling logger.setLevel(...) directly: the native core relies on its own "
		"level to decide whether to even format a message, for performance, so the two must "
		"stay in sync." );

	py::class_<EndoproteaseModel>(m, "EndoproteaseModel",
		"A single simulation: a generic endoprotease cutting one or more proteins, given "
		"per-position cleavage-probability data. See the project README for the JSON/dict "
		"configuration schema.")
		.def( py::init<>() )

		.def( "read_json", &EndoproteaseModel::readJson, py::arg("path"),
			"Load the simulation configuration from a JSON file on disk." )

		.def( "read_config", []( EndoproteaseModel& self, const nlohmann::json& config )
			{
				return self.readJsonObject( config );
			}, py::arg("config"),
			"Load the simulation configuration from a native Python dict (same schema as the "
			"JSON file format)." )

		.def( "run", &EndoproteaseModel::run, py::call_guard<py::gil_scoped_release>(),
			"Run the simulation to completion. Releases the GIL for the duration." )

		.def( "write_log", &EndoproteaseModel::writeLog, py::arg("path"),
			"Write the full simulation history to a CSV file (same format as the CLI's "
			"--output)." )

		.def( "compute_time_series", []( EndoproteaseModel& self, unsigned int period )
			{
				EndoproteaseModel::TimeSeries series = self.computeTimeSeries( period );

				py::dict result;
				result["time"] = series.time;
				result["time2"] = series.time2;
				result["enzyme"] = series.enzyme;
				for( size_t i = 0; i < series.peptideNames.size(); i++ )
					result[py::str(series.peptideNames[i])] = series.peptideCounts[i];

				return result;
			}, py::arg("period") = 10,
			"Compute the simulation's results in memory, as a flat dict of column name -> list "
			"of values (\"time\", \"time2\", \"enzyme\", plus one entry per distinct peptide "
			"produced during the simulation) -- ready to pass directly to pandas.DataFrame(...)." )

		// simulation parameters -- see the README's \"Configuration format\" section
		.def_readwrite( "max_time", &EndoproteaseModel::maxTime )
		.def_readwrite( "max_dh", &EndoproteaseModel::maxDH )
		.def_readwrite( "max_attempts", &EndoproteaseModel::maxAttempts )
		.def_readwrite( "max_attempts_per_time", &EndoproteaseModel::maxAttemptsPerTime )
		.def_readwrite( "current_enzyme", &EndoproteaseModel::currentEnzyme )
		.def_readwrite( "enzyme_always_dying", &EndoproteaseModel::enzymeAlwaysDying )
		.def_readwrite( "enzyme_dying_ratio", &EndoproteaseModel::enzymeDyingRatio )
		.def_readwrite( "random_seed", &EndoproteaseModel::randomSeed )
		.def_readonly( "t", &EndoproteaseModel::t,
			"Number of simulation iterations run so far." )
		;
}
