# Add the example application source file to a new executable target
add_executable(SplineCalcExample
    SplineCalcExample.cc
)

# Link the executable against our SplineCalculator library.
# CMake will automatically handle finding the library, its headers,
# and its dependencies (like ROOT).
target_link_libraries(SplineCalcExample PRIVATE SplineCalculator)