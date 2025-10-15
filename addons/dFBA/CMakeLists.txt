add_library(dFBA OBJECT)

file(GLOB SOURCES "src/*.cpp")
target_sources(dFBA PRIVATE ${SOURCES})
