set(BioInterface_SOURCES
    ${PROJECT_SOURCE_DIR}/bio_interface/Bio_microenvironment_interface.h
)

add_library(BioInterface INTERFACE)
target_sources(BioInterface INTERFACE ${BioInterface_SOURCES})
target_include_directories(BioInterface INTERFACE ${PROJECT_SOURCE_DIR})
