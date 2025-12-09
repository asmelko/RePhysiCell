get_filename_component(project ${CMAKE_CURRENT_SOURCE_DIR} NAME)
add_executable(${project})

file(GLOB MODULES "custom_modules/*.cpp")

target_sources(${project} PRIVATE main_ecoli_acetic_switch.cpp ${MODULES})

target_link_libraries(${project} PRIVATE BioFVM_utils BioFVM PhysiCell-core-no-addon
                                         PhysiCell-modules-no-addon)

target_include_directories(${project}
                           PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
target_include_directories(${project} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../..)

add_custom_command(
  TARGET ${project}
  POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory_if_different
          ${CMAKE_CURRENT_SOURCE_DIR}/config ${CMAKE_CURRENT_BINARY_DIR}/config
  COMMAND ${CMAKE_COMMAND} -E copy -t $<TARGET_FILE_DIR:${project}>
          $<TARGET_RUNTIME_DLLS:${project}>
  COMMAND_EXPAND_LISTS)

target_link_libraries(${project} PRIVATE dFBA)
