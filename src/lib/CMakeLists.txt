# PROJECT_DIR/SRC/LIB

cmake_minimum_required(VERSION 2.6)
project(SEM_LIB)

add_subdirectory(utilities)
add_subdirectory(iomanagment)
add_subdirectory(components)
add_subdirectory(time)
add_subdirectory(elements)
add_subdirectory(mesh)
add_subdirectory(fields)

add_subdirectory(solver)
add_subdirectory(fieldMath)
add_subdirectory(materials)
add_subdirectory(postprocessing)