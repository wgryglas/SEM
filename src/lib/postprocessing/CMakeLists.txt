# PROJECT_DIR/SRC/LIB/FIELDS
############# USER SETTINGS ####################################
set(TARGET postprocessing)

file(GLOB COMP_FILES *.cpp *.h)

set(LIB_NAMES fields fieldMath utilities)
set(THIRD_PARTY_NAMES boost)

############# BUILD #############################################
build_shared_lib(TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)


