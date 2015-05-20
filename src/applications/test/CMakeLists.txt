# PROJECT_DIR/SRC/APPLICATIONS/TEST
############## USER SETTINGS ####################################
set(TARGET SEMTest)

# set(COMP_FILES parseSEM.cpp)
# set(COMP_FILES generalTest.cpp)
# set(COMP_FILES writeField.cpp)
# set(COMP_FILES readField.cpp)
# set(COMP_FILES meshRead.cpp)
 set(COMP_FILES numArrays.cpp)
# set(COMP_FILES geometricField.cpp)
# set(COMP_FILES solveLaplace.cpp)
# set(COMP_FILES materialTest.cpp)
# set(COMP_FILES readingPath.cpp)
# set(COMP_FILES gmeshGeneratedMesh.cpp)
# set(COMP_FILES laplace_test.cpp)
# set(COMP_FILES vectorComponent.cpp)
# set(COMP_FILES div_grad.cpp)
# set(COMP_FILES elementFiledConversion.cpp)
# set(COMP_FILES etArray.cpp)
# set(COMP_FILES ode.cpp)
#set(COMP_FILES valueInterpolation.cpp)

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
