set(lib formfactor)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

file(GLOB src_files *.cpp)
set(api_files Polyhedron.h)

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
    FILES ${api_files}
    DESTINATION ${destination}/include
    COMPONENT Headers)
