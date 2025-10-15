add_library(BioFVM OBJECT)

file(GLOB SOURCES "*.cpp")
target_sources(BioFVM PRIVATE ${SOURCES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(BioFVM PUBLIC OpenMP::OpenMP_CXX)
else()
  message(
    WARNING "OpenMP not found. BioFVM will be built without OpenMP support.")
endif()
