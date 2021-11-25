set(lib formfactor)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

set(src_files)# formfactor.c)
set(inc_files)# formfactor.h)

add_library(${lib} ${src_files})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${formfactor_VERSION}
    SOVERSION ${formfactor_VERSION_MAJOR})

install(
    TARGETS ${lib}
    LIBRARY DESTINATION ${destination}/lib
    RUNTIME DESTINATION ${destination}/lib
    ARCHIVE DESTINATION ${destination}/lib
    COMPONENT Libraries)
install(
    FILES ${inc_files}
    DESTINATION ${destination}/include
    COMPONENT Headers)
