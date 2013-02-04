ADD_EXECUTABLE(rheometeranalyze main.cpp)
#===========================================================
TARGET_LINK_LIBRARIES(rheometeranalyze ${Boost_LIBRARIES})
INSTALL(TARGETS rheometeranalyze DESTINATION "bin/")
#===========================================================

