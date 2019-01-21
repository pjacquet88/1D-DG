cmake_minimum_required (VERSION 2.8)

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Green       "${Esc}[32m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
endif()

project(FWI)
enable_language(Fortran)

# Libraries
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Define source files
file(GLOB SRC "src/*.f90")

#--------------------- APPLICATIONS ---------------------------

# FORWARD APPLICATION
#set sources
set(FORWARD_SOURCES applications/forward.f90 ${SRC})
add_executable(FORWARD ${FORWARD_SOURCES})
#target compilation
target_compile_options(FORWARD PRIVATE -freal-4-real-8)
target_compile_options(FORWARD PRIVATE -fdefault-real-8)
#Add libraries
if(LAPACK_FOUND AND BLAS_FOUND)
  target_link_libraries(FORWARD ${BLAS_LIBRARIES})
  target_link_libraries(FORWARD ${LAPACK_LIBRARIES})
endif()
message(" ${BoldGreen}APPLICATION FORWARD CREATED !!! ${ColourReset}")


# FWI APPLICATION
#set sources
set(FWI_SOURCES applications/fwi.f90 ${SRC})
add_executable(FWI ${FWI_SOURCES})
#target compilation
target_compile_options(FWI PRIVATE -freal-4-real-8)
target_compile_options(FWI PRIVATE -fdefault-real-8)
#Add libraries
if(LAPACK_FOUND AND BLAS_FOUND)
  target_link_libraries(FWI ${BLAS_LIBRARIES})
  target_link_libraries(FWI ${LAPACK_LIBRARIES})
endif()
message(" ${BoldGreen}APPLICATION FWI CREATED !!! ${ColourReset}")


#--------------------- TESTS ---------------------------
enable_testing()

# TEST1 : Hello world
set(TEST1 tests/test_hello.f90)
add_executable(test1 ${TEST1})
add_test(Hello_word ${CMAKE_BINARY_DIR}/test1)

# TEST2 : Forward_without_param
add_test(Forward_without_param ${CMAKE_BINARY_DIR}/FORWARD)

# TEST3 : FWI_without_param
add_test(FWI_without_param_word ${CMAKE_BINARY_DIR}/FWI)
