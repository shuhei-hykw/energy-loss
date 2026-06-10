cmake_minimum_required(VERSION 3.16)
project(energy_loss_g4_table_generator CXX)

# Build configuration mirrors the install at ~/software/geant4/11.4.1
# (Multithreaded, GDML, Qt, OpenGL, RayTracer). Most of those are not
# strictly required for this minimal table generator, but matching the
# install avoids dynamic-loader surprises.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Single executable: reads CLI args (particle, material, energies, output)
# and writes a CSV with stopping power and CSDA range.
add_executable(g4_table_generator src/main.cc)
target_link_libraries(g4_table_generator PRIVATE ${Geant4_LIBRARIES})
