add_library(BioFVM OBJECT)
add_library(BioFVM_utils OBJECT)

target_sources(
  BioFVM
  PRIVATE BioFVM_agent_container.cpp
          BioFVM_basic_agent.cpp
          BioFVM_mesh.cpp
          BioFVM_microenvironment_adapter.cpp
          BioFVM_microenvironment.cpp
          BioFVM_MultiCellDS.cpp
          BioFVM_solvers.cpp)

target_sources(
  BioFVM_utils
  PRIVATE BioFVM_matlab.cpp
          BioFVM_utilities.cpp
          BioFVM_vector.cpp
          pugixml.cpp)

find_package(OpenMP 4)
if(OpenMP_CXX_FOUND)
  target_link_libraries(BioFVM PUBLIC OpenMP::OpenMP_CXX)
else()
  message(
    WARNING "OpenMP not found. BioFVM will be built without OpenMP support.")
endif()
