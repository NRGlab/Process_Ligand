cmake_minimum_required(VERSION 3.30)
project(Process_Ligand_real)

set(CMAKE_CXX_STANDARD 14)

# --- Static Build Option ---
# Enabled with -DBUILD_STATIC_EXECUTABLE=ON
option(BUILD_STATIC_EXECUTABLE "Build a fully static executable" OFF)

if(BUILD_STATIC_EXECUTABLE)
    message(STATUS "Static build enabled.")
    if(MSVC)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    elseif(UNIX AND NOT APPLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static")
    elseif(APPLE)
        message(WARNING "Full static linking (-static) is not supported on macOS and will not be enabled.")
    endif()
endif()

# --- User-configurable Open Babel path ---
# cmake -DSTATIC_INSTALL_DIR=/path/to/your/openbabel ..
set(STATIC_INSTALL_DIR "" CACHE PATH "Path to the static Open Babel installation directory")

if(NOT STATIC_INSTALL_DIR)
    message(FATAL_ERROR "STATIC_INSTALL_DIR is not set. Please specify the path to your Open Babel static installation using -DSTATIC_INSTALL_DIR=/path/to/openbabel")
endif()

include_directories(${STATIC_INSTALL_DIR}/include/openbabel3)

add_executable(Process_Ligand_real
        geometry.cpp
        Process_Ligand.cpp
        Process_Ligand.h)

target_link_directories(Process_Ligand_real PRIVATE ${STATIC_INSTALL_DIR}/bin)

if(WIN32)
    set(OPENBABEL_LIB_NAME openbabel-3)
else()
    set(OPENBABEL_LIB_NAME openbabel)
endif()

target_link_libraries(Process_Ligand_real PRIVATE ${OPENBABEL_LIB_NAME})
