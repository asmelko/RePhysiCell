add_library(PhysiCell-mechanics OBJECT)

target_sources(
  PhysiCell-mechanics
  PRIVATE PhysiCell_mechanics_agent.cpp
          PhysiCell_mechanics_standard_models.cpp
          PhysiCell_mechanics_environment.cpp
          PhysiCell_mechanics_functions.cpp)

find_package(OpenMP REQUIRED)
target_link_libraries(PhysiCell-mechanics PUBLIC OpenMP::OpenMP_CXX)
