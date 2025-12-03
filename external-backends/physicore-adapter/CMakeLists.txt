add_library(physicore_biofvm_adapter)

set(PHYSICORE_VERSION "1.2.0")

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  set(PHYSICORE_PACKAGE_URL
      https://github.com/bsc-life/PhysiCore/releases/download/v${PHYSICORE_VERSION}/physicore-linux-x64-v${PHYSICORE_VERSION}.tar.gz
  )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Darwin")
  set(PHYSICORE_PACKAGE_URL
      https://github.com/bsc-life/PhysiCore/releases/download/v${PHYSICORE_VERSION}/physicore-macos-arm64-v${PHYSICORE_VERSION}.tar.gz
  )
elseif(CMAKE_SYSTEM_NAME STREQUAL "Windows")
  set(PHYSICORE_PACKAGE_URL
      https://github.com/bsc-life/PhysiCore/releases/download/v${PHYSICORE_VERSION}/physicore-windows-x64-v${PHYSICORE_VERSION}.zip
  )
else()
  message(FATAL_ERROR "Unsupported OS: ${CMAKE_SYSTEM_NAME}")
endif()

include(FetchContent)
FetchContent_Declare(
  physicore
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  URL ${PHYSICORE_PACKAGE_URL})
FetchContent_MakeAvailable(physicore)

set(BioFVM_DIR ${physicore_SOURCE_DIR}/lib/cmake/BioFVM)
message(STATUS "Setting BioFVM_DIR to ${BioFVM_DIR}")

find_package(BioFVM CONFIG REQUIRED)

find_package(OpenMP 4)
if(OpenMP_CXX_FOUND)
  target_link_libraries(physicore_biofvm_adapter PUBLIC OpenMP::OpenMP_CXX)
else()
  message(
    WARNING "OpenMP not found. BioFVM will be built without OpenMP support.")
endif()

target_link_libraries(physicore_biofvm_adapter PUBLIC PhysiCore::reactions-diffusion.biofvm)

target_sources(physicore_biofvm_adapter PRIVATE
               basic_agent_adapter.cpp
               mesh_adapter.cpp
               microenvironment_adapter.cpp)

if(PHYSICELL_USE_PHYSICORE_BIOFVM)
  add_library(BioFVM ALIAS physicore_biofvm_adapter)
endif()
