cmake_minimum_required(VERSION 3.15)
project(PROTEIN_CUTTING CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# nlohmann/json is vendored as a single header under thirdparty/, so no
# network access or package manager is required to build (config parsing
# used to be done with tinyxml, also vendored in-tree, until the switch to JSON)
add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE ${PROTEIN_CUTTING_SOURCE_DIR}/thirdparty)

# spdlog is vendored header-only (just the include/ directory, no compiled lib),
# used for leveled logging instead of raw cout/cerr calls scattered through the code.
# SPDLOG_FMT_EXTERNAL is intentionally NOT defined, so spdlog uses its own bundled
# copy of fmt (spdlog/fmt/bundled/) rather than requiring a separately vendored fmt.
add_library(spdlog_headeronly INTERFACE)
target_include_directories(spdlog_headeronly INTERFACE ${PROTEIN_CUTTING_SOURCE_DIR}/thirdparty/spdlog)

# core simulation logic, as a library so it can be reused by both the CLI
# executable below and the Python bindings
add_library(protein-cutting-core STATIC EndoproteaseModel.cpp Peptide.cpp Log.cpp)
target_include_directories(protein-cutting-core PUBLIC ${PROTEIN_CUTTING_SOURCE_DIR})
target_link_libraries(protein-cutting-core PUBLIC nlohmann_json spdlog_headeronly)
# needed so this static library can be linked into the Python extension module (a shared
# library) below; a no-op on toolchains (e.g. MSVC) where it doesn't apply
set_target_properties(protein-cutting-core PROPERTIES POSITION_INDEPENDENT_CODE ON)

# CLI executable
add_executable(protein-cutting main.cpp)
target_link_libraries(protein-cutting PRIVATE protein-cutting-core)

# Python bindings (endocleave._core): optional, so the CLI above still builds standalone for
# anyone without a Python development environment. pybind11 itself is not vendored, unlike
# nlohmann/json and spdlog above -- it's a normal Python *build-time* dependency, resolved
# via "pip install pybind11" (or automatically, once packaged with scikit-build-core, from
# pyproject.toml's [build-system] requires).
find_package(Python COMPONENTS Interpreter Development QUIET)
find_package(pybind11 CONFIG QUIET)

if(Python_FOUND AND pybind11_FOUND)
	pybind11_add_module(_core bindings.cpp)
	target_link_libraries(_core PRIVATE protein-cutting-core)
	# pybind11_json is vendored header-only, same style as nlohmann/json and spdlog above
	target_include_directories(_core PRIVATE ${PROTEIN_CUTTING_SOURCE_DIR}/thirdparty/pybind11_json)

	if(MINGW)
		# statically link the MinGW runtime (libgcc/libstdc++/libwinpthread) into the module.
		# Without this, the module depends on those as external DLLs -- fine for the CLI
		# executable (built and run in the same MinGW shell where they're on PATH), but not
		# for a Python extension: since Python 3.8, the interpreter no longer searches PATH
		# for an extension module's DLL dependencies (a deliberate security hardening), so
		# importing the module would fail with a bare "DLL load failed" unless every consumer
		# manually called os.add_dll_directory(...). Static linking sidesteps this (and is
		# also what's needed for a redistributable wheel later).
		target_link_options(_core PRIVATE -static)
	endif()

	# scikit-build-core's convention: whatever is install()'d ends up in the wheel, rooted
	# at the install prefix it sets up; installing straight into "endocleave" merges the
	# compiled module into the same package directory as the pure-Python files it copies
	# from src/endocleave/ (see pyproject.toml's [tool.scikit-build] wheel.packages).
	install(TARGETS _core DESTINATION endocleave)
else()
	message(STATUS "Python and/or pybind11 not found: skipping the _core Python extension module (the CLI is unaffected)")
endif()
