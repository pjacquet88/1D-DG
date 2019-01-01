cmake_minimum_required (VERSION 2.8)
project (FWI)
enable_language(Fortran) 

#add the subdirectory
#add_subdirectory(lib)

# add the executable
set(FWI_SOURCES applications/fwi.f90 src/m_acoustic.f90 src/m_animation.f90 src/m_init_application.f90 src/m_powermethod.f90 src/m_adjoint.f90 src/m_file_function.f90 src/m_matrix.f90 src/m_time_scheme.f90 src/m_adjoint_test.f90 src/m_fwi.f90 src/m_polynom.f90)


add_executable(FWI ${FWI_SOURCES})

target_compile_options(FWI PRIVATE -freal-4-real-8)
target_compile_options(FWI PRIVATE -fdefault-real-8)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND AND BLAS_FOUND)
  message("BLAS LAPACK TROUVE !!!")
  target_link_libraries(FWI ${BLAS_LIBRARIES})
  target_link_libraries(FWI ${LAPACK_LIBRARIES})
endif()
