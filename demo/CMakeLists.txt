cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(demo_using_formfactor VERSION 0.1.2 LANGUAGES CXX)

## Compiler settings.
option(WERROR "Treat warnings as errors" OFF)

set(CMAKE_CXX_STANDARD 17)
if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    message(STATUS "compiling for Windows")
else()
    option(PEDANTIC "Compile with pedantic warnings" ON)
    if(PEDANTIC)
        add_compile_options(-pedantic -Wall)
    endif()
endif()
add_compile_options(-O2)
if(WERROR)
    add_compile_options(-Werror)
endif()

## Dependencies
find_package(LibHeinz REQUIRED CONFIG)
find_package(formfactor REQUIRED CONFIG)

## Source
set(demos
    octahedron
    )

## Configure executable
foreach(app ${demos})
    add_executable(${app} "${app}.cpp")
    target_link_libraries(${app} formfactor)
endforeach()
