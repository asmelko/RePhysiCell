add_library(PhysiCell-modules OBJECT)
add_library(PhysiCell-modules-no-addon OBJECT)

set(PHYSICELL_MODULES_ADDON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/PhysiCell_MultiCellDS.cpp)

target_sources(
  PhysiCell-modules
  PRIVATE PhysiCell_geometry.cpp
          PhysiCell_MultiCellDS.cpp
          PhysiCell_pathology.cpp
          PhysiCell_POV.cpp
          PhysiCell_pugixml.cpp
          PhysiCell_settings.cpp
          PhysiCell_SVG.cpp
          PhysiCell_various_outputs.cpp)

target_sources(
  PhysiCell-modules-no-addon
  PRIVATE PhysiCell_geometry.cpp
          PhysiCell_pathology.cpp
          PhysiCell_POV.cpp
          PhysiCell_pugixml.cpp
          PhysiCell_settings.cpp
          PhysiCell_SVG.cpp
          PhysiCell_various_outputs.cpp)
