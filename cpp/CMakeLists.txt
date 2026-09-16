cmake_minimum_required(VERSION 3.15)
project(PROTEIN_CUTTING CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
	# spdlog bundles fmt, which refuses to compile under MSVC at all without this: a
	# static_assert enforces the choice explicitly rather than silently assuming the wrong
	# source/execution character set. Found by CI (this project has no MSVC build locally
	# to catch it against) -- GCC/Clang default to UTF-8 and never needed this flag.
	add_compile_options(/utf-8)
endif()

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

if(MINGW)
	# statically link the MinGW runtime into the CLI too, for the same reason as the Python
	# module below: a dynamically-linked binary depends on libstdc++/libgcc/libwinpthread as
	# external DLLs, resolved via the OS's normal DLL search order at startup -- which is a
	# real problem the moment more than one MinGW toolchain exists on the same machine (e.g.
	# Git for Windows bundles its own MinGW64, alongside a separately installed one): if a
	# *different*, ABI-incompatible libstdc++-6.dll is found first, the result is a silent,
	# hard-to-diagnose crash inside the C++ standard library at runtime, in code that
	# compiled and linked without any error. Static linking removes the dependency (and
	# hence the ambiguity) entirely.
	target_link_options(protein-cutting PRIVATE -static)
endif()

# Python bindings (endocleave._core): optional, so the CLI above still builds standalone for
# anyone without a Python development environment. pybind11 itself is not vendored, unlike
# nlohmann/json and spdlog above -- it's a normal Python *build-time* dependency, resolved
# via "pip install pybind11" (or automatically, once packaged with scikit-build-core, from
# pyproject.toml's [build-system] requires).
#
# NOTE: deliberately NOT also calling find_package(Python ...) ourselves first.
# pybind11_add_module() already does its own Python detection internally -- and does it more
# reliably than a separate, explicit find_package(Python COMPONENTS Interpreter Development)
# call turned out to be: inside cibuildwheel's manylinux containers specifically, that separate
# call silently failed to find the container's per-version Python, even though it's the exact
# same Python pybind11 itself found a moment later on its own. Since the build was gated on
# *both* being found, the whole block was silently skipped -- producing a "successful" wheel
# with no compiled extension inside it at all (caught by cibuildwheel's test-command, and by
# auditwheel refusing to repair a wheel with no shared library in it). Gating on pybind11_FOUND
# alone avoids re-implementing (and re-breaking) what pybind11 already gets right.
find_package(pybind11 CONFIG QUIET)

if(pybind11_FOUND)
	pybind11_add_module(_core bindings.cpp)
	target_link_libraries(_core PRIVATE protein-cutting-core)
	# pybind11_json is vendored header-only, same style as nlohmann/json and spdlog above
	target_include_directories(_core PRIVATE ${PROTEIN_CUTTING_SOURCE_DIR}/thirdparty/pybind11_json)

	if(MINGW)
		# see the CLI target above for the general reason (avoiding a possible mismatched
		# libstdc++-6.dll at runtime); on top of that, a Python extension module has its own,
		# stricter problem: since Python 3.8, the interpreter no longer searches PATH for an
		# extension module's DLL dependencies at all (a deliberate security hardening), so
		# importing the module would fail outright with a bare "DLL load failed" unless every
		# consumer manually called os.add_dll_directory(...). Static linking sidesteps both
		# issues at once (and is also what's needed for a redistributable wheel later).
		target_link_options(_core PRIVATE -static)
	endif()

	# scikit-build-core's convention: whatever is install()'d ends up in the wheel, rooted
	# at the install prefix it sets up; installing straight into "endocleave" merges the
	# compiled module into the same package directory as the pure-Python files it copies
	# from src/endocleave/ (see pyproject.toml's [tool.scikit-build] wheel.packages).
	install(TARGETS _core DESTINATION endocleave)
else()
	message(STATUS "pybind11 not found: skipping the _core Python extension module (the CLI is unaffected)")
endif()
