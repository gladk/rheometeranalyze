ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp particle.cpp rheometer.cpp band.cpp export.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES} vtkHybrid)
INSTALL(TARGETS rheometeranalyze DESTINATION "bin/")
#===========================================================

