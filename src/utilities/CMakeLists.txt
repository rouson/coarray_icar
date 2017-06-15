set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

#Toggle C preprocessor macro for turning assertions on or off
if(NO_ASSERTIONS)
  set_source_files_properties(assertions_interface.f90 PROPERTIES COMPILE_FLAGS "-cpp -DUSE_ASSERTIONS=.false.")
else()
  set_source_files_properties(assertions_interface.f90 PROPERTIES COMPILE_FLAGS "-cpp -DUSE_ASSERTIONS=.true.")
endif()

add_library( icar_utils
  assertions_implementation.f90
  assertions_interface.f90
  configuration_implementation.f90
  configuration_interface.f90
)

install(TARGETS icar_utils DESTINATION lib)

# Find the native NetCDF includes and library
include(${CMAKE_SOURCE_DIR}/cmake-scripts/cmake-modules/FindNetCDF.cmake)

set (NETCDF_F90 "YES") # Require Fortran interfaces to be FOUND
find_package (NetCDF REQUIRED)

include_directories(${NETCDF_INCLUDES})

# For fortran module files it's a good practice to use CMake's "object libraries"
# to avoid parallel build race conditions. See
# https://cmake.org/cmake/help/latest/command/add_library.html#object-libraries

add_library(model_const OBJECT from-GARD/model_constants.f90)

add_library(io_routines from-GARD/io_routines.f90 $<TARGET_OBJECTS:model_const>)

target_link_libraries (io_routines ${NETCDF_LIBRARIES})
