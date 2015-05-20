# PROJECT_DIR/SRC/LIB/MATEIRALS
############# USER SETTINGS ####################################
set(TARGET materials)

file(GLOB COMP_FILES *.cpp *.h)

set(LIB_NAMES)
set(THIRD_PARTY_NAMES)

############# BUILD #############################################
build_shared_lib(TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)


