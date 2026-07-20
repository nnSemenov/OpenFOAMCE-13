
add_subdirectory(specie)
add_subdirectory(thermophysicalProperties)

add_subdirectory(basic)
add_subdirectory(multicomponentThermo)
add_subdirectory(chemistryModel)

add_subdirectory(solidThermo)

set(FOAM_regular_libraries ${FOAM_regular_libraries} PARENT_SCOPE)