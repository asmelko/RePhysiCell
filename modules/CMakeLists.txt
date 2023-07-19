cmake_minimum_required(VERSION 3.18)

project(PhysiCellModules VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Target PhysiCellModules
file(GLOB_RECURSE src_files "*")

add_library(PhysiCellModules OBJECT ${src_files})

if(MSVC)
    target_compile_options(PhysiCellModules PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(PhysiCellModules PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()
