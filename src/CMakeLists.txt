ADD_EXECUTABLE(rheometeranalyze main.cpp config.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES})
INSTALL(TARGETS rheometeranalyze DESTINATION "bin/")
#===========================================================

