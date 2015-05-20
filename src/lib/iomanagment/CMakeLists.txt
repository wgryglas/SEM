# PROJECT_DIR/SRC/LIB/IOMANAGMENT
############# USER SETTINGS ####################################
#------------ NAME ---------------------------------------------
#set(TARGET iomanagment)

#------------ INCLUDE PATHS -------------------------------------
#set(PUBLIC_INCLUDE)# ${LIB_INCLUDE})
#set(PRIVATE_INCLUDE)# ${LIB_SRC}/iomanagment )
#set(INCLUDE
#        ${PUBLIC_INCLUDE}
#        ${PRIVATE_INCLUDE}
#        ${THIRD_PARTY_INCLUDE} )

#------------ LIBRARIES -----------------------------------------
#set(LIB boost_filesystem)
#set(LIB_PATH ${THIRD_PARTY_BIN}./boost)

#------------ FILES ---------------------------------------------
#file(GLOB SOURCE_FILES *.cpp *.h)
#file(GLOB INCLUDE_FILES ${PUBLIC_INCLUDE}/iomanagment/*h ${PRIVATE_INCLUDE}/*h)
#set(FILES ${SOURCE_FILES} ${INCLUDE_FILES} )
#set(FILES ${SOURCE_FILES})


#build_shared_lib(TARGET FILES INCLUDE LIB LIB_PATHS)

# PROJECT_DIR/SRC/LIB/IOMANAGMENT
############# USER SETTINGS ####################################
set(TARGET iomanagment)

file(GLOB COMP_FILES *.cpp *.h)

set(LIB_NAMES
       boost_filesystem
       boost_system
       boost_regex)

set(THIRD_PARTY_NAMES boost)

############# BUILD #############################################
build_shared_lib(TARGET COMP_FILES LIB_NAMES THIRD_PARTY_NAMES)


