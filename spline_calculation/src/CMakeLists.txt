# Add the SplineCalculator source files to a new library target
add_library(SplineCalculator SHARED
    SplineCalculator.cxx
    SplineCalculator.h
)

# Link the library against the required ROOT libraries.
# ROOT::DataFrame will pull in all necessary dependencies.
target_link_libraries(SplineCalculator PUBLIC ${ROOT_LIBRARIES})

# Include ROOT headers
target_include_directories(SplineCalculator PUBLIC 
    ${ROOT_INCLUDE_DIRS}
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)