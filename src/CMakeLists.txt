ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp particle.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES})
INSTALL(TARGETS rheometeranalyze DESTINATION "bin/")
#===========================================================

