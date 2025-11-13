cmake_minimum_required(VERSION 3.16)
project(SplineCalculation)

# Set a standard C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the ROOT package, which is a core dependency for the library.
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RooFit RooFitCore MathCore RIO Net) # Include necessary components

# Include ROOT's use file to set up proper compilation flags and include directories
include(${ROOT_USE_FILE})

# Add the subdirectory that builds our library
add_subdirectory(src)

# Add the subdirectory that builds our example application
add_subdirectory(app)