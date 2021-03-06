project(RheometerAnalyze C CXX)
cmake_minimum_required(VERSION 3.18)

INCLUDE(FindOpenMP)
INCLUDE(GNUInstallDirs)
#===========================================================
IF(NOT(CMAKE_CXX_FLAGS))
  #If flags are not set, add all usefull flags
  SET(CMAKE_BUILD_TYPE Release)
  SET(CMAKE_CXX_FLAGS  "-Wall -fPIC -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
  set (CMAKE_SHARED_LINKER_FLAGS "-Wl,--as-needed")
ENDIF (NOT(CMAKE_CXX_FLAGS))

IF(OPENMP_FOUND)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DENABLED_OPENMP ${OpenMP_CXX_FLAGS}")
  MESSAGE(STATUS "Found OpenMP")
ELSE(OPENMP_FOUND)
  MESSAGE(STATUS "OpenMP NOT found")
ENDIF(OPENMP_FOUND)

#C+11 support
set(CMAKE_CXX_FLAGS "-std=c++0x ${CMAKE_CXX_FLAGS}")
#===========================================================
# Add possibility to use local boost installation (e.g. -DLocalBoost=1.46.1)

SET(BOOST_MIN_VERSION "1.54.0")
FIND_PACKAGE(Boost ${LocalBoost} COMPONENTS program_options filesystem system regex iostreams log log_setup thread REQUIRED)
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
ADD_DEFINITIONS( "-DBOOST_LOG_DYN_LINK" )

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})

FIND_PACKAGE(VTK COMPONENTS CommonCore IOImage IOXML REQUIRED)
FIND_PACKAGE(ALGLIB)
#===========================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")
#===========================================================

