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
#  Reader.cxx Flux.cxx FluxManager.cxx Calculator.cxx Writer.cxx DirWriter.cxx SRWriter.cxx
  Reader.cxx SampleManager.cxx ObservableManager.cxx
  )

set(includes_IO
#  Reader.h Flux.h FluxManager.h Calculator.h Writer.h DirWriter.h SRWriter.h Data.h Constants.h
  Reader.h SampleManager.h ObservableManager.h
  )

############################Weightor############################

add_library(Utils_IO SHARED ${source_IO} ${includes_IO})

target_include_directories(Utils_IO PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
  "$<INSTALL_INTERFACE:include/Utils_IO>")

install(TARGETS Utils_IO
  EXPORT Utils_IOTargets
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES DESTINATION include/Utils_IO)



ROOT_GENERATE_DICTIONARY(stl_dict MODULE Utils_IO LINKDEF LinkDef.h)

# Install the header files.
install(FILES ${includes_IO} DESTINATION include/Utils_IO)

# Make sure the current directories are available for the later
# compilation.  This is required to make ROOT_GENERATE_DICTIONARY
# work.
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# target_include_directories(Utils_IO PUBLIC OscLib::All Eigen3::Eigen)
target_link_libraries(Utils_IO LINK_PUBLIC OscLib::All Eigen3::Eigen
                      duneanaobj::StandardRecordProxy
                      duneanaobj::StandardRecordFlat
                      duneanaobj::StandardRecord
                      duneanaobj::StandardRecord_dict
		      yaml-cpp)
target_include_directories(Utils_IO
  PRIVATE
  $ENV{SRPROXY_INC}
)
