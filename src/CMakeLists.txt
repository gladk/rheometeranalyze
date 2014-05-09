#===========================================================
ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp particle.cpp rheometer.cpp 
  band.cpp export.cpp force.cpp bandRow.cpp forceRow.cpp snapshot.cpp 
  snapshotRow.cpp particleRow.cpp math_custom.cpp interactionori.cpp
  forceChain.cpp)
  
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES})
INSTALL(TARGETS rheometeranalyze DESTINATION ${CMAKE_INSTALL_BINDIR})
#===========================================================
IF(${VTK_MAJOR_VERSION} EQUAL 6)
  TARGET_LINK_LIBRARIES(rheometeranalyze ${VTK_LIBRARIES})
  MESSAGE(STATUS "VTK version 6 is found")
  ADD_DEFINITIONS("-DVTK6")
ELSE(${VTK_MAJOR_VERSION} EQUAL 6)
  TARGET_LINK_LIBRARIES(rheometeranalyze vtkHybrid)
ENDIF(${VTK_MAJOR_VERSION} EQUAL 6)
#===========================================================
IF(ALGLIB_FOUND)
  ADD_DEFINITIONS(-DALGLIB)
  INCLUDE_DIRECTORIES(${ALGLIB_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES(rheometeranalyze ${ALGLIB_LIB})
  MESSAGE(STATUS "ALGLIB is found and shearband export will be produced")
ELSE(ALGLIB_FOUND)
  MESSAGE(STATUS "ALGLIB is NOT found and shearband export will NOT be produced")
ENDIF(ALGLIB_FOUND)
#===========================================================

