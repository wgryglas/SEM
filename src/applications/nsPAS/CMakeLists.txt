# PROJECT_DIR/SRC/APPLICATIONS/LAPLACE
############## USER SETTINGS ####################################
set(TARGET nsPAS)

set(COMP_FILES nsPAS.cpp)

set(LIB_NAMES
    boost_system
    boost_filesystem
    boost_regex
    iomanagment
    components
    time
    elements
    mesh
    las
    fields
    fieldMath
    materials
    utilities)

set(THIRD_PARTY_NAMES boost)

############## BUILD #############################################
build_application(TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)
