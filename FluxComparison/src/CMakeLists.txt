message("TPCSim Library")

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release)
# endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS "-Werror -Wall -Wpedantic")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

###################################################

# Configure the dependencies
find_package(ROOT REQUIRED COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  cmessage(STATUS "ROOT_USEFILE: ${ROOT_USE_FILE}")
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Define the source and include files that should be used for the io
# part of ERepSim.

set(source   FluxReader.cxx BartolFluxReader.cxx HondaFluxReader.cxx)
set(includes FluxReader.h   BartolFluxReader.h   HondaFluxReader.h  )

############################FluxComparison############################

add_library(FluxComparison SHARED ${source} ${includes})

target_link_libraries(FluxComparison yaml-cpp)

target_include_directories(FluxComparison PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/FluxComparison>")

install(TARGETS FluxComparison
  EXPORT FluxComparisonIOTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/FluxComparison)

ROOT_GENERATE_DICTIONARY(stl_dict MODULE FluxComparison LINKDEF LinkDef.h)

# Install the header files.
install(FILES ${includes} DESTINATION include/FluxComparison)

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(FluxComparison
  PRIVATE
  $ENV{SRPROXY_INC}
)
