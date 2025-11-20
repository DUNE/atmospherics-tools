cmake_minimum_required(VERSION 3.14)

# Set the policy version minimum to allow older dependencies (like yaml-cpp) to configure.
# This is the equivalent of passing -DCMAKE_POLICY_VERSION_MINIMUM=3.5 on the command line.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

project(SplineCalculation)

# Set a default build type if none is specified.
# Common options: Debug, Release, RelWithDebInfo, MinSizeRel
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the build type" FORCE)
  message(STATUS "No build type specified, defaulting to ${CMAKE_BUILD_TYPE}")
endif()

# Set a standard C++ version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(DUNE_ANAOBJ_BRANCH v03_11_00)

# Find the ROOT package, which is a core dependency for the library.
find_package(ROOT REQUIRED COMPONENTS Core Hist Tree RooFit RooFitCore MathCore RIO Net) # Include necessary components

# Include ROOT's use file to set up proper compilation flags and include directories
include(${ROOT_USE_FILE})

#=================== CPM    

# download CPM.cmakes
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.2/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

# Add yaml-cpp using CPM. It will first try find_package, then download if not found.
CPMAddPackage(
  NAME yaml-cpp
  VERSION 0.8.0
  GIT_TAG 0.8.0
  GITHUB_REPOSITORY jbeder/yaml-cpp
  GIT_SHALLOW YES
  OPTIONS
    "YAML_CPP_INSTALL ON"
    "YAML_CPP_BUILD_TESTS OFF"
    "YAML_CPP_BUILD_CONTRIB OFF"
    "YAML_BUILD_SHARED_LIBS ON"
)
if(NOT TARGET yaml-cpp::yaml-cpp)
  cmessage(FATAL_ERROR "yaml-cpp package not found and could not be installed.")
endif()

find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO EG)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
else()
  message(FATAL_ERROR "ROOT not found!")
endif(ROOT_FOUND)

CPMAddPackage(
 	  NAME duneanaobj
	  GIT_TAG ${DUNE_ANAOBJ_BRANCH}
	  GITHUB_REPOSITORY DUNE/duneanaobj 
	  DOWNLOAD_ONLY YES
  )
	include_directories(${duneanaobj_SOURCE_DIR})

  ROOT_GENERATE_DICTIONARY(StandardRecordDict ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/StandardRecord.h ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/SRGlobal.h
  LINKDEF ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/classes_def.xml)

  file(GLOB SR_IMPL_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.cxx")
  LIST(APPEND SR_IMPL_FILES ${CMAKE_CURRENT_BINARY_DIR}/StandardRecordDict.cxx)

  add_library(duneanaobj_StandardRecord SHARED ${SR_IMPL_FILES})
  target_link_libraries(duneanaobj_StandardRecord PUBLIC ${ROOT_LIBRARIES})

  target_include_directories(duneanaobj_StandardRecord PUBLIC 
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}> 
    $<INSTALL_INTERFACE:include>)
  target_include_directories(duneanaobj_StandardRecord PRIVATE 
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord> #root puts this in the dictionary
    )
  
  set_target_properties(duneanaobj_StandardRecord PROPERTIES EXPORT_NAME all)
  install(TARGETS duneanaobj_StandardRecord EXPORT mach3dune-targets DESTINATION lib)
  
  file(GLOB SR_HEADER_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.h")

  install(FILES ${SR_HEADER_FILES} DESTINATION include/duneanaobj/StandardRecord)
  install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict_rdict.pcm 
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict.rootmap 
    DESTINATION lib)

  add_library(duneanaobj::all ALIAS duneanaobj_StandardRecord)

# Add the subdirectory that builds our library
add_subdirectory(src)

# Add the subdirectory that builds our example application
add_subdirectory(app)