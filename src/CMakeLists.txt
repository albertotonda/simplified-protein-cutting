cmake_minimum_required(VERSION 3.15)
project(PROTEIN_CUTTING CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# nlohmann/json is vendored as a single header under thirdparty/, so no
# network access or package manager is required to build (config parsing
# used to be done with tinyxml, also vendored in-tree, until the switch to JSON)
add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE ${PROTEIN_CUTTING_SOURCE_DIR}/thirdparty)

# core simulation logic, as a library so it can be reused by both the CLI
# executable below and the future Python bindings
add_library(protein-cutting-core STATIC PepsinModel.cpp Peptide.cpp)
target_include_directories(protein-cutting-core PUBLIC ${PROTEIN_CUTTING_SOURCE_DIR})
target_link_libraries(protein-cutting-core PUBLIC nlohmann_json)

# CLI executable
add_executable(protein-cutting main.cpp)
target_link_libraries(protein-cutting PRIVATE protein-cutting-core)
