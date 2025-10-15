add_library(PhysiMeSS OBJECT)

file(GLOB SOURCES "*.cpp")
target_sources(PhysiMeSS PRIVATE ${SOURCES})
