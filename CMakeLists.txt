cmake_minimum_required (VERSION 2.8)
project(FWI)
enable_language(Fortran) 

# Library
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Define source file
file(GLOB SRC "src/*.f90")

# FWI APPLICATION
#set sources
set(FWI_SOURCES applications/fwi.f90 ${SRC})
#add_executable(FWI applications/fwi.f90 ${SRC})
add_executable(FWI ${FWI_SOURCES})
#target compilation
target_compile_options(FWI PRIVATE -freal-4-real-8)
target_compile_options(FWI PRIVATE -fdefault-real-8)
#Add libraries
if(LAPACK_FOUND AND BLAS_FOUND)
  message("BLAS LAPACK TROUVE !!!")
  target_link_libraries(FWI ${BLAS_LIBRARIES})
  target_link_libraries(FWI ${LAPACK_LIBRARIES})
endif()
