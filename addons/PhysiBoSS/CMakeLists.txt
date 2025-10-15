add_library(PhysiBoSS OBJECT)

target_sources(
  PhysiBoSS
  PRIVATE src/maboss_intracellular.cpp src/maboss_network.cpp
          ${PROJECT_SOURCE_DIR}/core/PhysiCell_cell.cpp
          ${PROJECT_SOURCE_DIR}/modules/PhysiCell_MultiCellDS.cpp)

set(MABOSS_VERSION v2.5.6)
set(MABOSS_MAX_NODES
    64
    CACHE STRING "Maximum number of nodes for MaBoSS")

if(NOT MABOSS_MAX_NODES EQUAL 64)
  set(MABOSS_LIB_NAME MaBoSS_${MABOSS_MAX_NODES}n-static)
else()
  set(MABOSS_LIB_NAME MaBoSS-static)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(MABOSS_ARCHIVE libMaBoSS-linux64.tar.gz)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_OSX_ARCHITECTURE STREQUAL
                                               "x86_64")
  set(MABOSS_ARCHIVE libMaBoSS-osx64.tar.gz)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin" AND CMAKE_OSX_ARCHITECTURE STREQUAL
                                               "arm64")
  set(MABOSS_ARCHIVE libMaBoSS-macos-arm64.tar.gz)
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(MABOSS_ARCHIVE libMaBoSS-win64.zip)
else()
  message(FATAL_ERROR "Unsupported OS: ${CMAKE_SYSTEM_NAME}")
endif()

include(FetchContent)
FetchContent_Declare(
  maboss
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL https://github.com/sysbio-curie/MaBoSS/releases/download/${MABOSS_VERSION}/${MABOSS_ARCHIVE}
)
FetchContent_MakeAvailable(maboss)

target_include_directories(PhysiBoSS PUBLIC ${maboss_SOURCE_DIR}/engine/include)
target_link_directories(PhysiBoSS PUBLIC ${maboss_SOURCE_DIR}/engine/lib)
target_link_libraries(PhysiBoSS PUBLIC ${MABOSS_LIB_NAME})
target_compile_definitions(PhysiBoSS PUBLIC MAXNODES=${MABOSS_MAX_NODES}
                                            ADDON_PHYSIBOSS)
