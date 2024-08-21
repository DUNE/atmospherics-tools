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
find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
cmessage(STATUS "ROOT_USEFILE: ${ROOT_USE_FILE}")
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

find_package(OscLib REQUIRED $ENV{OSCLIB_INC} $ENV{OSCLIB_LIB})
find_package(Eigen3Local REQUIRED $ENV{EIGEN_INC})
cmessage(STATUS "EIGEN: ${OscLib_LIBS}")

# Define the source and include files that should be used for the io
# part of ERepSim.

set(source_IO
  Reader.cxx Flux.cxx FluxManager.cxx Calculator.cxx Writer.cxx DirWriter.cxx
  )

set(includes_IO
  Reader.h Flux.h FluxManager.h Calculator.h Writer.h DirWriter.h
  )

############################Weightor############################

add_library(WeightorIO SHARED ${source_IO} ${includes_IO})

target_include_directories(WeightorIO PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/WeightorIO>")

install(TARGETS WeightorIO
  EXPORT WeightorIOTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/WeightorIO)



ROOT_GENERATE_DICTIONARY(stl_dict MODULE WeightorIO LINKDEF LinkDef.h)

# Install the header files.
install(FILES ${includes_IO} DESTINATION include/WeightorIO)

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# target_include_directories(WeightorIO PUBLIC OscLib::All Eigen3::Eigen)
target_link_libraries(WeightorIO LINK_PUBLIC OscLib::All Eigen3::Eigen
                      duneanaobj::StandardRecordProxy
                      duneanaobj::StandardRecordFlat
                      duneanaobj::StandardRecord
                      duneanaobj::StandardRecord_dict)
target_include_directories(WeightorIO
  PRIVATE
  $ENV{SRPROXY_INC}
)
