cmake_minimum_required(VERSION 2.6)
project(PROTEIN_CUTTING)

add_subdirectory(tinyxml)

# Set include directories
include_directories(${PROTEIN_CUTTING_SOURCE_DIR})
include_directories(${PROTEIN_CUTTING_SOURCE_DIR}/tinyxml)

# Add executable
add_executable(protein-cutting main.cpp PepsinModel.cpp Peptide.cpp)

# Link the executable with the freshly created library
target_link_libraries(protein-cutting Tinyxml)
