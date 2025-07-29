cmake_minimum_required(VERSION 3.15)
project(ProcessLigand LANGUAGES CXX)

# --- Static Build Option ---
# Enabled by passing -DBUILD_STATIC_EXECUTABLE=ON to the cmake command 
option(BUILD_STATIC_EXECUTABLE "Build a fully static executable" OFF)

if(BUILD_STATIC_EXECUTABLE)
    message(STATUS "Static build enabled.")
    if(MSVC)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
endif()

# --- User-configurable Open Babel Path ---
# Allow the user to specify the path to Open Babel when running cmake (path to folder containing bin, lib, include),...:
# cmake -DSTATIC_INSTALL_DIR=/path/to/your/openbabel .. 
set(STATIC_INSTALL_DIR "" CACHE PATH "Path to the static Open Babel installation directory")

if(NOT STATIC_INSTALL_DIR)
    message(FATAL_ERROR "STATIC_INSTALL_DIR is not set. Please specify the path using -DSTATIC_INSTALL_DIR=/path/to/openbabel")
endif()

add_executable(ProcessLigand
        geometry.cpp
        Process_Ligand.cpp
)

set_target_properties(ProcessLigand PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED ON
)

target_include_directories(ProcessLigand PRIVATE
        ${STATIC_INSTALL_DIR}/include/openbabel3
)

if(WIN32)
    set(OPENBABEL_LIB_DIR ${STATIC_INSTALL_DIR}/bin)
    set(OPENBABEL_LIB_NAME openbabel-3)
else()
    set(OPENBABEL_LIB_DIR ${STATIC_INSTALL_DIR}/lib)
    set(OPENBABEL_LIB_NAME openbabel)
endif()

target_link_directories(ProcessLigand PRIVATE ${OPENBABEL_LIB_DIR})
target_link_libraries(ProcessLigand PRIVATE ${OPENBABEL_LIB_NAME})

if(BUILD_STATIC_EXECUTABLE)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(ProcessLigand PRIVATE -static)
    elseif(APPLE)
        message(WARNING "Full static linking (-static) is not supported on macOS and will be skipped.")
    endif()
endif()
