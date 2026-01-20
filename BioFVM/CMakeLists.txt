add_library(BioFVM OBJECT)
add_library(BioFVM_utils OBJECT)

target_sources(
  BioFVM
  PRIVATE BioFVM_agent_container.cpp
          BioFVM_basic_agent.cpp
          BioFVM_microenvironment.cpp
          BioFVM_solvers.cpp
          BioFVM_implementation.cpp
          BioFVM_legacy_implementation.cpp)

target_sources(
  BioFVM_utils
  PRIVATE BioFVM_vector.cpp
          BioFVM_mesh.cpp
          BioFVM_basic_agent_PIMPL.cpp
          BioFVM_utilities.cpp
          BioFVM_matlab.cpp
          BioFVM_MultiCellDS.cpp
          pugixml.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(BioFVM PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(BioFVM_utils PUBLIC OpenMP::OpenMP_CXX)
