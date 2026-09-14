// Thin wrapper around a single, shared, named spdlog logger, used throughout the
// protein-cutting core instead of raw cout/cerr calls.
//
// The core library only ever logs through this; it never decides *where* messages
// end up (console, file, or later -- for the Python bindings -- a sink that forwards
// records into Python's own `logging` module). That choice belongs to whoever links
// the core: see log::setSinks() / log::setLevel(), which main.cpp calls based on the
// CLI flags.

#ifndef __LOG__
#define __LOG__

#include <spdlog/spdlog.h>
#include <sstream>
#include <vector>
#include <memory>

// NOTE: this is named "logging", not "log" -- most of this codebase has a blanket
// "using namespace std;", and std::log(double) (the natural logarithm, from <cmath>)
// would otherwise shadow a namespace called "log" and break "log::get()" et al.
namespace logging
{
	// name of the shared logger; also the name a future Python `logging.getLogger()` would use
	constexpr const char* LOGGER_NAME = "protein_cutting";

	// the shared logger instance; silent (null sink) until setSinks() is called
	std::shared_ptr<spdlog::logger>& get();

	// replace the logger's sinks entirely (e.g. console, or console + file)
	void setSinks( std::vector<spdlog::sink_ptr> sinks );

	// set the minimum level that will actually be logged/formatted
	void setLevel( spdlog::level::level_enum level );
}

// Stream-style logging macros, kept close to the "cout << a << b << c" style already
// used throughout this codebase, but routed through spdlog. The should_log() check
// means the message is only ever formatted (the ostringstream work) if that level is
// actually enabled, so a disabled LOG_TRACE(...) costs almost nothing at runtime --
// this replaces both the old "verbose" bool and the hardcoded "DEBUG" macro.
#define LOG_AT(lvl, expr) \
	do { \
		if( logging::get()->should_log(lvl) ) { \
			std::ostringstream _log_stream; \
			_log_stream << expr; \
			logging::get()->log(lvl, _log_stream.str()); \
		} \
	} while(0)

#define LOG_TRACE(expr) LOG_AT(spdlog::level::trace, expr)
#define LOG_DEBUG(expr) LOG_AT(spdlog::level::debug, expr)
#define LOG_INFO(expr)  LOG_AT(spdlog::level::info, expr)
#define LOG_WARN(expr)  LOG_AT(spdlog::level::warn, expr)
#define LOG_ERROR(expr) LOG_AT(spdlog::level::err, expr)

#endif
