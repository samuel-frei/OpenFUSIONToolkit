install(FILES oft_gs.in DESTINATION examples/MUG/equilibrium COMPONENT examples)
configure_file( oft_gs.in ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
oft_add_exe(driver_taylor.F90 )
