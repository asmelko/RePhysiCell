cmake_minimum_required(VERSION 3.18)

project(BioFVM VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP REQUIRED)

# Target BioFVM
file(GLOB_RECURSE src_files "*")

add_library(BioFVM OBJECT ${src_files})

if(MSVC)
    target_compile_options(BioFVM PRIVATE $<$<COMPILE_LANGUAGE:CXX>:/bigobj>)
else()
    target_compile_options(BioFVM PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-march=native>)
endif()
