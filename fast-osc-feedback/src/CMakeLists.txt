add_library(Oscillogram SHARED Oscillogram.cxx)

message(STATUS "Include dir -> ${OSCPROB_INCLUDE_DIR}")
target_include_directories(Oscillogram PUBLIC ${OSCPROB_INCLUDE_DIR})
ROOT_GENERATE_DICTIONARY(Oscillogram_dict Oscillogram.h MODULE Oscillogram LINKDEF LinkDef.h)
target_link_libraries(Oscillogram PRIVATE ROOT::RIO)