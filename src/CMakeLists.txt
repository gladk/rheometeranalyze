IF(ALGLIB_FOUND)
  ADD_DEFINITIONS(-DALGLIB)
ENDIF(ALGLIB_FOUND)

#===========================================================
ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp particle.cpp rheometer.cpp band.cpp export.cpp force.cpp bandRow.cpp forceRow.cpp snapshot.cpp snapshotRow.cpp particleRow.cpp math_custom.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES} ${ALGLIB_LIB} vtkHybrid)

INSTALL(TARGETS rheometeranalyze DESTINATION ${CMAKE_INSTALL_BINDIR})
#===========================================================

