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

add_executable(weightor weightor.cxx)
target_link_libraries(weightor
  LINK_PUBLIC WeightorIO ROOT::EG ${ROOT_LIBRARIES} OscLib argparse)
target_compile_features(weightor PRIVATE cxx_std_17)

target_include_directories(weightor
  PRIVATE
  $ENV{SRPROXY_INC}
)

add_executable(split_channels split_channels.cxx)
target_link_libraries(split_channels
  LINK_PUBLIC WeightorIO ROOT::EG ${ROOT_LIBRARIES})
target_include_directories(split_channels
  PRIVATE
  $ENV{SRPROXY_INC}
)

install(TARGETS weightor split_channels RUNTIME DESTINATION bin)

