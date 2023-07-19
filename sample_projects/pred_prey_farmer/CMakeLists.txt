cmake_minimum_required(VERSION 3.18)

project(pred_prey_farmer VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(pred_prey_farmer main.cpp custom_modules/custom.cpp)
target_include_directories(pred_prey_farmer PRIVATE "../..")
target_link_libraries(pred_prey_farmer PRIVATE BioFVM PhysiCellCore PhysiCellModules OpenMP::OpenMP_CXX)

if(MSVC)
    target_compile_options(pred_prey_farmer PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(pred_prey_farmer PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

add_custom_command(TARGET pred_prey_farmer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/config $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)

add_custom_command(TARGET pred_prey_farmer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:${PROJECT_NAME}>/output)