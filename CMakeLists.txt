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

# Define source files
file(GLOB SRC "src/*90")

# Libraries
add_library(fwi_lib ${SRC})
message("Library created")
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
target_compile_options(fwi_lib PRIVATE -DKIND_REAL=sp)

#--------------------- APPLICATIONS ---------------------------

# FORWARD APPLICATION
#set sources
set(FORWARD_SOURCES applications/forward.f90)
add_executable(FORWARD ${FORWARD_SOURCES})
#target compilation

#Add libraries
target_link_libraries(FORWARD fwi_lib)
if(LAPACK_FOUND AND BLAS_FOUND)
  target_link_libraries(FORWARD ${BLAS_LIBRARIES})
  target_link_libraries(FORWARD ${LAPACK_LIBRARIES})
endif()
message(" ${BoldGreen}APPLICATION FORWARD CREATED !!! ${ColourReset}")


#FWI APPLICATION
#set sources
set(FWI_SOURCES applications/fwi.f90)
add_executable(FWI ${FWI_SOURCES})
#target compilation

#Add libraries
target_link_libraries(FWI fwi_lib)
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


#------------------ CUSTOM TARGETS --------------------
add_custom_target(data_clean
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/Files/*.dat
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/build/fort.*
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/build/*.gif
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/build/*.dat
  COMMAND ${CMAKE_COMMAND} -E remove ${CMAKE_SOURCE_DIR}/build/*.a)

add_custom_target(cleanall
  COMMAND make clean
  COMMAND make data_clean)
