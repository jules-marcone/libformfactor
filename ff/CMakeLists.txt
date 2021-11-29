set(lib formfactor)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

file(GLOB src_files *.cpp)
set(api_files Polyhedron.h PolyhedralTopology.h PolyhedralComponents.h)

add_library(${lib} ${src_files})

target_include_directories(${lib}
    PRIVATE
    ${LibHeinz_INCLUDE_DIR}
    PUBLIC
    ${CMAKE_SOURCE_DIR})

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
    COMPONENT Libraries)

install(
    FILES ${api_files}
    DESTINATION include/ff
    COMPONENT Headers)
