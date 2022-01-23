set(demos
    tetrahedron
    )

foreach(app ${demos})
    add_executable(${app} ${app}.cpp)
#    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_link_libraries(${app} formfactor)
endforeach()
