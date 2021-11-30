file(GLOB test_sources "*Test.cpp")

if (WIN32)
    add_compile_options(/wd4068 /wd4623 /wd4625 /wd4626 /wd4710 /wd4711 /wd4820 )
endif()

foreach(source_file ${test_sources})
    get_filename_component(base_name ${source_file} NAME_WE)
    add_executable(${base_name} ${source_file} catch2.cpp)
    target_include_directories(${base_name}
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR} ${LibHeinz_INCLUDE_DIR})
    target_link_libraries(${base_name} ${formfactor_LIBRARY})
    add_test(NAME ${base_name} COMMAND ${base_name})
endforeach(source_file)
