add_library(libRoadrunner OBJECT)

file(GLOB SOURCES "src/*.cpp")
target_sources(libRoadrunner PRIVATE ${SOURCES})
