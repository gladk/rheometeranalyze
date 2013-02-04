project(RheometerAnalyze C CXX)
cmake_minimum_required(VERSION 2.8)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cMake")

INCLUDE(FindOpenMP)
#===========================================================
IF(NOT(CMAKE_CXX_FLAGS))
  #If flags are not set, add all usefull flags
  SET(CMAKE_BUILD_TYPE Release)
  SET(CMAKE_CXX_FLAGS  "-Wall -fPIC -g -O2 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Wformat-security -Werror=format-security")
ENDIF (NOT(CMAKE_CXX_FLAGS))
#===========================================================
# Add possibility to use local boost installation (e.g. -DLocalBoost=1.46.1)

FIND_PACKAGE(Boost ${LocalBoost} COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES (${Boost_INCLUDE_DIRS})

FIND_PACKAGE(Eigen3 REQUIRED)
INCLUDE_DIRECTORIES(${EIGEN2_INCLUDE_DIR})

#===========================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/src")
#===========================================================

