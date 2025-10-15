add_library(PhysiCell-core OBJECT)
add_library(PhysiCell-core-no-addon OBJECT)

set(PHYSICELL_CORE_ADDON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/PhysiCell_cell.cpp)

target_sources(
  PhysiCell-core
  PRIVATE PhysiCell_basic_signaling.cpp
          PhysiCell_cell_container.cpp
          PhysiCell_cell.cpp
          PhysiCell_constants.cpp
          PhysiCell_custom.cpp
          PhysiCell_phenotype.cpp
          PhysiCell_rules.cpp
          PhysiCell_signal_behavior.cpp
          PhysiCell_standard_models.cpp
          PhysiCell_utilities.cpp)

target_sources(
  PhysiCell-core-no-addon
  PRIVATE PhysiCell_basic_signaling.cpp
          PhysiCell_cell_container.cpp
          PhysiCell_constants.cpp
          PhysiCell_custom.cpp
          PhysiCell_phenotype.cpp
          PhysiCell_rules.cpp
          PhysiCell_signal_behavior.cpp
          PhysiCell_standard_models.cpp
          PhysiCell_utilities.cpp)
