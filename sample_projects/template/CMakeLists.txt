cmake_minimum_required(VERSION 3.18)

project(template VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

add_executable(template main.cpp custom_modules/custom.cpp)
target_include_directories(template PRIVATE "../..")
target_link_libraries(template PRIVATE BioFVM PhysiCellCore PhysiCellModules OpenMP::OpenMP_CXX)

if(MSVC)
    target_compile_options(template PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(template PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()

add_custom_command(TARGET template POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    ${CMAKE_CURRENT_SOURCE_DIR}/config $<TARGET_FILE_DIR:${PROJECT_NAME}>/config)

add_custom_command(TARGET template POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory
    $<TARGET_FILE_DIR:${PROJECT_NAME}>/output)