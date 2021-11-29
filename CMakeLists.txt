cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(formfactor VERSION 0.1.0 LANGUAGES CXX)

## Options.

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
option(WERROR "Treat warnings as errors" OFF)

## Compiler settings.

set(CMAKE_CXX_STANDARD 17)
if(WIN32)
    set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
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

## Dependences.

find_package(LibHeinz REQUIRED)
message(STATUS "LibHeinz: found=${LibHeinz_FOUND}, include_dirs=${LibHeinz_INCLUDE_DIR}, "
    "version=${LibHeinz_VERSION}")

## Install settings.

# install(EXPORT FormFactorTargets FILE FormFactorTargets.cmake DESTINATION cmake)

## Subdirectories.

include(CTest)

add_subdirectory(ff)
add_subdirectory(test)

## CPack settings.

include(InstallRequiredSystemLibraries)
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_PACKAGE_VERSION_MAJOR "${formfactor_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${formfactor_VERSION_MINOR}")
set(CPACK_SOURCE_GENERATOR "TGZ")
include(CPack)
include(CMakePackageConfigHelpers)
export(PACKAGE formfactor)

## Config files.

configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/formfactorConfig.cmake"
    INSTALL_DESTINATION "lib/cmake/example"
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/formfactorConfigVersion.cmake"
    VERSION "${formfactor_VERSION_MAJOR}.${formfactor_VERSION_MINOR}"
    COMPATIBILITY AnyNewerVersion
    )

install(FILES
    ${PROJECT_BINARY_DIR}/formfactorConfig.cmake
    ${PROJECT_BINARY_DIR}/formfactorConfigVersion.cmake
    DESTINATION cmake COMPONENT dev)
