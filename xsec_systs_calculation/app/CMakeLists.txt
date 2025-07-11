LIST(APPEND TARGETS_TO_BUILD 
UpdateReweight
ReadCAF
PrintWeights
)

foreach(targ ${TARGETS_TO_BUILD})
  add_executable(${targ} ${targ}.cxx)
  target_link_libraries(${targ} 
    nusyst::systproviders
    ${ROOT_LIBRARIES}
    ROOT::Geom
    ROOT::MathMore
    duneanaobj::StandardRecordProxy
    duneanaobj::StandardRecordFlat
    duneanaobj::StandardRecord
    duneanaobj::StandardRecord_dict
  )
  target_include_directories(${targ}
    PRIVATE
    $ENV{SRPROXY_INC}
  )

endforeach()

install(TARGETS ${TARGETS_TO_BUILD} DESTINATION bin/)
