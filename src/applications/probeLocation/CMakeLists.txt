# PROJECT_DIR/SRC/APPLICATIONS/probeLocation
############## USER SETTINGS ####################################
set(TARGET probeLocation)

set(COMP_FILES probeLocation.cpp)

set(LIB_NAMES
    boost_system
    boost_filesystem
    boost_regex
    iomanagment
    time
    materials
    components
    elements
    mesh
    las
    fields
    fieldMath
    utilities
    postprocessing)

set(THIRD_PARTY_NAMES boost)

############## BUILD #############################################
build_application(TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)
