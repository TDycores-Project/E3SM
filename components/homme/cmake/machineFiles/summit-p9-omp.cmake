#interactive job
#bsub -W 2:00 -nnodes 1 -P cli115 -Is /bin/bash

#cmake -C ~/acme-fork-lb/components/homme/cmake/machineFiles/summit.cmake -DHOMMEXX_MPI_ON_DEVICE=FALSE ~/acme-fork-lb/components/homme/
#SET (HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")

SET (NETCDF_DIR $ENV{OLCF_NETCDF_FORTRAN_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{OLCF_HDF5_ROOT} CACHE FILEPATH "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET(WITH_PNETCDF FALSE CACHE FILEPATH "")

SET(USE_QUEUING FALSE CACHE BOOL "")

SET(ENABLE_CUDA FALSE CACHE BOOL "")

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET(BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET(HOMME_ENABLE_COMPOSE FALSE CACHE BOOL "")
SET(HOMMEXX_GB_CONFIG ON CACHE BOOL "")

SET (HOMMEXX_VECTOR_SIZE 8 CACHE STRING "")
SET (AVX_VERSION 0 CACHE STRING "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA OFF CACHE BOOL "")
SET(Kokkos_ARCH_POWER9 ON CACHE BOOL "")
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "mpicc" CACHE STRING "")

SET (ADD_Fortran_FLAGS "-fopenmp -O3" CACHE STRING "")
SET (ADD_C_FLAGS       "-fopenmp -O3 --std=c++11" CACHE STRING "")
SET (ADD_CXX_FLAGS     "-fopenmp -O3 --std=c++11" CACHE STRING "")

SET(LINKER_ADD_ON "-lstdc++" CACHE STRING "")

set (ENABLE_OPENMP ON CACHE BOOL "")
set (ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
set (ENABLE_HORIZ_OPENMP ON CACHE BOOL "")

set (HOMME_TESTING_PROFILE "dev" CACHE STRING "")

set (USE_NUM_PROCS 4 CACHE STRING "")

set (OPT_FLAGS " -mcpu=power9 -mtune=power9 -DNDEBUG " CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to core" CACHE FILEPATH "")
