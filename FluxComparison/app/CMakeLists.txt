# Build all of the applications here.
project(app)

set(CMAKE_CXX_FLAGS "-Werror -Wall -Wpedantic")

# Configure the dependencies                                                                                                                                                                                          
find_package(ROOT REQUIRED COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  cmessage(STATUS "ROOT_USEFILE: ${ROOT_USE_FILE}")
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

add_executable(FluxComp FluxComp.cxx)
target_link_libraries(FluxComp LINK_PUBLIC FluxComparison ROOT::Tree ROOT::Hist ROOT::Physics ${ROOT_LIBRARIES})
target_compile_features(FluxComp PRIVATE cxx_std_17)

target_include_directories(FluxComp
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS FluxComp RUNTIME DESTINATION bin)

