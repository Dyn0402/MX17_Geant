cmake_minimum_required(VERSION 3.16 FATAL_ERROR)
project(MicromegasSim VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Geant4 -- all physics modules
find_package(Geant4 REQUIRED ui_all vis_all)

# Use Geant4 include dirs and definitions
include(${Geant4_USE_FILE})

# Find ROOT (optional but recommended for output)
find_package(ROOT QUIET)

# Include dirs
include_directories(${PROJECT_SOURCE_DIR}/include)

# Source files
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cc)

# Build the executable
add_executable(mm_sim mm_sim.cc ${SOURCES})
target_link_libraries(mm_sim ${Geant4_LIBRARIES})

if(ROOT_FOUND)
    target_include_directories(mm_sim PRIVATE ${ROOT_INCLUDE_DIRS})
    target_link_libraries(mm_sim ${ROOT_LIBRARIES})
    target_compile_definitions(mm_sim PRIVATE USE_ROOT)
    message(STATUS "ROOT found and enabled: ${ROOT_VERSION}")
else()
    message(STATUS "ROOT not found; CSV output only")
endif()

# Copy macros to build dir
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})

install(TARGETS mm_sim DESTINATION bin)
