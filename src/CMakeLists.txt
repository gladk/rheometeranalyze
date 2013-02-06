ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp particle.cpp rheometer.cpp band.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES})
INSTALL(TARGETS rheometeranalyze DESTINATION "bin/")
#===========================================================

