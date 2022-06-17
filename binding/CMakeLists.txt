cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

project(binding_attempt)

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

if(WERROR)
    add_compile_options(-Werror)
endif()


find_package(PythonLibs)
include_directories(${PYTHON_INCLUDE_DIRS})
link_libraries(${PYTHON_LIBRARIES})
#link_libraries(${PYTHON_DEBUG_LIBRARIES})

include_directories(${PYTHON_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

find_package(LibHeinz REQUIRED CONFIG)
find_package(formfactor REQUIRED CONFIG)

find_package(SWIG REQUIRED)
include(${SWIG_USE_FILE})

set(CMAKE_SWIG_OUTDIR ${CMAKE_CURRENT_BINARY_DIR}/python)
set(SWIG_OUTFILE_DIR ${CMAKE_CURRENT_BINARY_DIR}/cpp)

set_property(SOURCE pave.i PROPERTY CPLUSPLUS ON)

swig_add_library(pave TYPE SHARED LANGUAGE python SOURCES pave.i pave.cpp)

swig_link_libraries(pave ${PYTHON_LIBRARIES})