set(lib formfactor)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

file(GLOB src_files *.cpp)
set(api_files Polyhedron.h Prism.h PolyhedralTopology.h PolyhedralComponents.h)

add_library(${lib} ${src_files})

target_include_directories(${lib}
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
    )
target_include_directories(${lib} PRIVATE ${LibHeinz_INCLUDE_DIR})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${formfactor_VERSION}
    SOVERSION ${formfactor_VERSION_MAJOR})

install(
    TARGETS ${lib}
    EXPORT formfactorTargets
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION lib
    ARCHIVE DESTINATION lib)

install(
    FILES ${api_files}
    DESTINATION include/ff)
