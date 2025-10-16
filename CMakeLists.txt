cmake_minimum_required(VERSION 3.12)
project(
  PhysiCell
  LANGUAGES CXX
  VERSION 1.14.2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

add_subdirectory(BioFVM)
add_subdirectory(core)
add_subdirectory(custom_modules)
add_subdirectory(modules)
add_subdirectory(examples)

add_subdirectory(addons/dFBA)
add_subdirectory(addons/libRoadrunner)
add_subdirectory(addons/PhysiBoSS)
add_subdirectory(addons/PhysiMeSS)

add_subdirectory(sample_projects/asymmetric_division)
add_subdirectory(sample_projects/biorobots)
add_subdirectory(sample_projects/cancer_biorobots)
add_subdirectory(sample_projects/celltypes3)
add_subdirectory(sample_projects/custom_division)
add_subdirectory(sample_projects/episode)
add_subdirectory(sample_projects/heterogeneity)
add_subdirectory(sample_projects/immune_function)
add_subdirectory(sample_projects/interactions)
add_subdirectory(sample_projects/mechano)
add_subdirectory(sample_projects/physimess)
add_subdirectory(sample_projects/pred_prey_farmer)
add_subdirectory(sample_projects/rules_sample)
add_subdirectory(sample_projects/template)
add_subdirectory(sample_projects/virus_macrophage)
add_subdirectory(sample_projects/worm)

add_subdirectory(sample_projects_intracellular/boolean/cancer_invasion)
add_subdirectory(sample_projects_intracellular/boolean/physiboss_cell_lines)
add_subdirectory(sample_projects_intracellular/boolean/template_BM)
add_subdirectory(sample_projects_intracellular/boolean/tutorial)

add_subdirectory(sample_projects_intracellular/fba/cancer_metabolism)
# broken target - can not compile
# add_subdirectory(sample_projects_intracellular/fba/ecoli_acetic_switch) 

add_subdirectory(sample_projects_intracellular/ode/ode_energy)
