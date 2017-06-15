include_directories(${CMAKE_BINARY_DIR}/mod)

add_executable( initialization-test
  test-initialization.f90
  ../objects/domain_implementation.f90
  ../objects/domain_interface.f90
  ../objects/exchangeable_implementation.f90
  ../objects/exchangeable_interface.f90
)
add_executable( plane-wave
  plane-wave.f90
  ../objects/domain_implementation.f90
  ../objects/domain_interface.f90
  ../objects/exchangeable_implementation.f90
  ../objects/exchangeable_interface.f90
)
target_link_libraries(initialization-test icar_utils)
target_link_libraries(plane-wave icar_utils io_routines)

configure_file(input-parameters.txt ${CMAKE_BINARY_DIR}/input-parameters.txt COPYONLY)

install(TARGETS initialization-test DESTINATION bin)
install(TARGETS plane-wave DESTINATION bin)
