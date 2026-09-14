// Implementation of the shared logger wrapper, see Log.h

#include <Log.h>
#include <spdlog/sinks/null_sink.h>

namespace logging
{
	std::shared_ptr<spdlog::logger>& get()
	{
		static std::shared_ptr<spdlog::logger> instance = []
		{
			// default to a null sink, so nothing is printed at all until a real sink
			// is configured (main() does this for the CLI; the Python bindings will
			// do the equivalent, forwarding records into Python's logging module)
			auto nullSink = std::make_shared<spdlog::sinks::null_sink_mt>();
			auto logger = std::make_shared<spdlog::logger>( LOGGER_NAME, nullSink );
			logger->set_level( spdlog::level::info );
			spdlog::register_logger( logger );
			return logger;
		}();

		return instance;
	}

	void setSinks( std::vector<spdlog::sink_ptr> sinks )
	{
		get()->sinks() = std::move(sinks);
	}

	void setLevel( spdlog::level::level_enum level )
	{
		get()->set_level( level );
	}
}
