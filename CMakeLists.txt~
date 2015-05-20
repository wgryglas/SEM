#PROJECT_DIR

cmake_minimum_required(VERSION 2.6)
project(SEM)

#################################################
#       BUILD SETUP
#################################################
# comand line argument for debug:
#-DCMAKE_BUILD_TYPE=Debug
#SET(CMAKE_CXX_FLAGS "-std=c++0x")
SET(CMAKE_CXX_FLAGS "-std=c++0x -O3")

#################################################
#       BUILD VARIABLES SETUP
#################################################
# ---------- Setup lib build output Directories -------------------------
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/out/lib )
# --------- Setup the build Executable output Directory -------------
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/out/application )

#################################################
#       GLOBAL VARIABLES
#################################################
set(PRJ /home/wojtek/cpp/SEM/SEM_PROJECT)
set(SYS_DIR_NAME linux)

set(SRC ${PRJ}/src)
set(THIRD_PARTY ${PRJ}/third-party)
#set(EXPORT ${PRJ}/export/${SYS_DIR_NAME})
set(EXPORT /home/wojtek/Apps/SMESH/SEM)

set(APP_SRC ${SRC}/application)
set(LIB_SRC ${SRC}/lib)
set(LIB_INCLUDE ${SRC}/include)
set(ETC ${SRC}/etc)
set(THIRD_PARTY_INCLUDE ${THIRD_PARTY}/include)

set(BUILD_APP_BIN ${CMAKE_BINARY_DIR}/out/application)
set(BUILD_LIB_BIN ${CMAKE_BINARY_DIR}/out/lib)

set(THIRD_PARTY_BIN ${THIRD_PARTY}/bin/${SYS_DIR_NAME})

set(EXPORT_APP_BIN ${EXPORT}/application)
set(EXPORT_LIB_BIN ${EXPORT}/lib)
set(EXPORT_INCLUDE ${EXPORT}/include)

#################################################
#       MACROS
#################################################
macro(build_application TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)

    cmake_minimum_required(VERSION 2.6)

    project(SEM_${${TARGET}}_APPLICATION)

    #include dirs
        # headers from libraries
    include_directories(${LIB_INCLUDE})
        #third party location
    include_directories(${THIRD_PARTY_INCLUDE})

    #libraries locations
        #compiled libraries
    message(${BUILD_LIB_BIN})
    link_directories(${BUILD_LIB_BIN})
        #third party compiled libraries
    foreach(IN ${${THIRD_PARTY_NAMES}})
        link_directories(${THIRD_PARTY_BIN}/${IN})
    endforeach(IN)

    add_executable(${${TARGET}} ${${COMP_FILES}})

    target_link_libraries(${${TARGET}} ${${LIB_NAMES}})

    install(TARGETS ${${TARGET}}
            RUNTIME DESTINATION ${EXPORT_APP_BIN})

endmacro(build_application)

macro(build_shared_lib TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)

    cmake_minimum_required(VERSION 2.6)

    project(SEM_${${TARGET}}_SHARED_LIB)

    #find relative location to "src/lib/..."
    file(RELATIVE_PATH REL_LOC ${LIB_SRC} ${CMAKE_CURRENT_SOURCE_DIR})

    #include dirs
        # current library public headers location-to allow including in current library
        # files without full local path to header, but just by header name
    include_directories(${LIB_INCLUDE}/${REL_LOC})
        # headers from internal project libraries
    include_directories(${LIB_INCLUDE})
        #third party headers location
    include_directories(${THIRD_PARTY_INCLUDE})

    #libraries locations
        #compiled libraries
    link_directories(${BUILD_LIB_BIN})
        #third party compiled libraries
    foreach(IN ${${THIRD_PARTY_NAMES}})
        link_directories(${THIRD_PARTY_BIN}/${IN})
    endforeach(IN)

    #build shared lib
    # add extra public lib headers - it will alow qtCreator to show in project those headers
    # if you would like to see also private headers, then you need to add them insied COMP_FILES variable
    file(GLOB PUBLIC_HEADERS ${LIB_INCLUDE}/${REL_LOC}/*.h)
    add_library(${${TARGET}} SHARED ${${COMP_FILES}} ${PUBLIC_HEADERS})

    #link libraries
    target_link_libraries(${${TARGET}} ${${LIB_NAMES}}) #link chosen libraries

    #install nesseccery inlude files
    install(FILES ${PUBLIC_HEADERS} DESTINATION ${EXPORT_INCLUDE}/${REL_LOC})
    #install compiled shared library
    install(TARGETS ${${TARGET}}
            RUNTIME DESTINATION ${EXPORT_LIB_BIN}
            LIBRARY DESTINATION ${EXPORT_LIB_BIN}
            ARCHIVE DESTINATION ${EXPORT_LIB_BIN})

endmacro(build_shared_lib)

macro(general_build_application TARGET FILES INCLUDE LIB LIB_PATHS)

    cmake_minimum_required(VERSION 2.6)

    project(SEM_${${TARGET}})

    include_directories(${${INCLUDE}})

    link_directories(${${LIB_PATHS}})

    add_executable(${${TARGET}} ${${FILES}})

    target_link_libraries(${${TARGET}} ${${LIB}})

endmacro(general_build_application)


macro(general_build_shared_lib TARGET FILES INCLUDE LIB LIB_PATHS)

    cmake_minimum_required(VERSION 2.6)

    project(SEM_${${TARGET}})

    include_directories(${${INCLUDE}})

    link_directories(${${LIB_PATHS}})

    add_library(${${TARGET}} SHARED ${${FILES}})

    target_link_libraries(${${TARGET}} ${${LIB}})

endmacro(general_build_shared_lib)


##################################################
# COMPILE SUBDIRECTORIES
##################################################
add_subdirectory(src)


##################################################
# EXPORT EXTRA DIRECTORIES
##################################################
install(DIRECTORY ${ETC} DESTINATION ${EXPORT})

