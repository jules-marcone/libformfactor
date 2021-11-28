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
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib
    COMPONENT Libraries
    EXPORT FormFactorTargets)
install(
    EXPORT FormFactorTargets
    FILE FormFactorTargets.cmake
    DESTINATION cmake)
install(
    FILES ${api_files}
    DESTINATION include
    COMPONENT Headers)
