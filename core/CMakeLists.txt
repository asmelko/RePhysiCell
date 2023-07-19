cmake_minimum_required(VERSION 3.18)

project(PhysiCellCore VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

# Target PhysiCellCore
file(GLOB_RECURSE src_files "*")
list(FILTER src_files EXCLUDE REGEX "PhysiCell_digital_cell_line.cpp")

add_library(PhysiCellCore OBJECT ${src_files})

if(MSVC)
    target_compile_options(PhysiCellCore PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(PhysiCellCore PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native -fomit-frame-pointer -mfpmath=both -m64 -fopenmp>)
endif()
