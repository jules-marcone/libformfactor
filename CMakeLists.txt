cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(libformfactor VERSION 0.1 LANGUAGES CXX)

## Options.

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)

## Compiler settings.

set(CMAKE_CXX_STANDARD 17)
if(WIN32)
    set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
    message(STATUS "compiling for Windows")
endif()
add_compile_options(-O2)
if(PEDANTIC)
    add_compile_options(-pedantic -Wall)
endif()
if(WERROR)
    add_compile_options(-Werror)
endif()

## Dependences.

find_package(LibHeinz REQUIRED)
message(STATUS "LibHeinz: found=${LibHeinz_FOUND}, include_dirs=${LibHeinz_INCLUDE_DIR}, "
    "version=${LibHeinz_VERSION}")

## Install settings.

configure_file(libformfactor.pc.in libformfactor.pc @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/libformfactor.pc DESTINATION lib/pkgconfig/)

include(CTest)

add_subdirectory(lib)
add_subdirectory(test)
