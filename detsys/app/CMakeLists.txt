# Build all of the applications here.
project(app)

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release)
# endif()

# set(CMAKE_CXX_FLAGS "-Wall -Wextra -O3 -funroll-loops -march=native")
set(CMAKE_CXX_FLAGS "-Werror -Wall -Wpedantic")

# Set the default verbosity of the compilation...
set(CMAKE_VERBOSE_MAKEFILE TRUE)

find_package(ROOT REQUIRED
  COMPONENTS Geom Physics Matrix MathCore Tree RIO)
if(ROOT_FOUND)
  include(${ROOT_USE_FILE})
endif(ROOT_FOUND)

add_executable(ComparisonScript ComparisonScript.cxx)
target_link_libraries(ComparisonScript
  LINK_PUBLIC Utils_IO ROOT::EG ${ROOT_LIBRARIES} GENIE3::All NuOscillator yaml-cpp)
target_compile_features(ComparisonScript PRIVATE cxx_std_17)

target_include_directories(ComparisonScript
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS ComparisonScript RUNTIME DESTINATION bin)

