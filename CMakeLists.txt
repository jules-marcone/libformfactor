cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(formfactor VERSION 0.1.2 LANGUAGES CXX)

## Options.

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
option(WERROR "Treat warnings as errors" OFF)

## Compiler settings.

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

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

## Dependences.

find_package(LibHeinz REQUIRED)
message(STATUS "LibHeinz: found=${LibHeinz_FOUND}, include_dirs=${LibHeinz_INCLUDE_DIR}, "
    "version=${LibHeinz_VERSION}")

## Subdirectories.

include(CTest)

add_subdirectory(ff)
add_subdirectory(test)
add_subdirectory(demo)

## CPack settings.

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${formfactor_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${formfactor_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)

## Install settings.

install(EXPORT formfactorTargets
  FILE formfactorTargets.cmake
  DESTINATION cmake)

## Config files.

include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/formfactorConfig.cmake"
    INSTALL_DESTINATION cmake
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO    )

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/formfactorConfigVersion.cmake"
    VERSION "${formfactor_VERSION_MAJOR}.${formfactor_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion    )

install(FILES
    ${PROJECT_BINARY_DIR}/formfactorConfig.cmake
    ${PROJECT_BINARY_DIR}/formfactorConfigVersion.cmake
    DESTINATION cmake)

## Generate export targets.

export(EXPORT formfactorTargets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/formfactorTargets.cmake)
