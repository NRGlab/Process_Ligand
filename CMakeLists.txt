cmake_minimum_required(VERSION 3.15)
project(ProcessLigand LANGUAGES CXX)

include(FetchContent)

option(BUILD_STATIC_EXECUTABLE "Build a statically linked executable" OFF)

if(BUILD_STATIC_EXECUTABLE)
    message(STATUS "Building a static executable")
    if(MSVC)
        set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    endif()
    set(BUILD_SHARED OFF CACHE BOOL "Force disable shared builds")
else()
    message(STATUS "Building a dynamic executable")
    set(BUILD_SHARED On CACHE BOOL "Force disable shared builds")
endif()

set(BUILD_GUI OFF CACHE BOOL "Disable building the Open Babel GUI")
set(PYTHON_BINDINGS OFF CACHE BOOL "Disable Python language bindings")
set(PERL_BINDINGS OFF CACHE BOOL "Disable Perl language bindings")
set(RUBY_BINDINGS OFF CACHE BOOL "Disable Ruby language bindings")
set(CSHARP_BINDINGS OFF CACHE BOOL "Disable C# language bindings")
set(JAVA_BINDINGS OFF CACHE BOOL "Disable Java language bindings")
set(WITH_MAEPARSER OFF CACHE BOOL "Disable the MAEParser feature")
set(PTHREAD_LIBRARY "" CACHE FILEPATH "Disable PThread library")
set(ZLIB_LIBRARY OFF CACHE FILEPATH "Disable ZLIB library")
set(ZLIB_INCLUDE_DIR OFF CACHE PATH "Disable ZLIB includes")
set(INCHI_LIBRARY OFF CACHE FILEPATH "Disable InChI library")
set(LIBXML2_LIBRARY OFF CACHE FILEPATH "Disable LibXML2 library")
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)

FetchContent_Declare(
        openbabel
        GIT_REPOSITORY https://github.com/ThomasDesc/openbabel.git
        GIT_TAG master
)

FetchContent_MakeAvailable(openbabel)

message(STATUS "Open Babel source directory: ${openbabel_SOURCE_DIR}")
message(STATUS "Open Babel binary directory: ${openbabel_BINARY_DIR}")

add_executable(ProcessLigand
        geometry.cpp
        Process_Ligand.cpp
)
set_property(TARGET ProcessLigand PROPERTY
        MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
target_link_libraries(ProcessLigand PRIVATE openbabel)

target_include_directories(ProcessLigand PRIVATE
        ${openbabel_SOURCE_DIR}/include
        ${openbabel_BINARY_DIR}/include
)

if(BUILD_STATIC_EXECUTABLE)
    if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        target_link_options(ProcessLigand PRIVATE -static)
    elseif(MINGW)
        target_link_options(ProcessLigand PRIVATE -static)
    elseif(APPLE)
        message(WARNING "Full static linking (-static) is not supported on macOS and will be skipped.")
    endif()
endif()