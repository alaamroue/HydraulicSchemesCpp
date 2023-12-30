cmake_minimum_required(VERSION 3.16)


set(CMAKE_CONFIGURATION_TYPES Debug;Release)

project(godunovCpp)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

set(CPP_H_FILES
    source_code/1_CLDomainCartesian.cpp
    source_code/1_CLDomainCartesian.h
    source_code/2_CLFriction.cpp
    source_code/2_CLFriction.h
    source_code/3_CLSolverHLLC.cpp
    source_code/3_CLSolverHLLC.h
    source_code/4_CLDynamicTimestep.cpp
    source_code/4_CLDynamicTimestep.h
    source_code/5_CLSchemeGodunov.cpp
    source_code/5_CLSchemeGodunov.h
    source_code/6_CLBoundaries.cpp
    source_code/6_CLBoundaries.h
    source_code/7_CLSchemePromaides.cpp
    source_code/7_CLSchemePromaides.h
    source_code/definitions.h
    source_code/GlobalHandlerClass.cpp
    source_code/GlobalHandlerClass.h
    source_code/globals_handlers.cpp
    source_code/main.cpp
    source_code/main.h
    source_code/normalPlain.cpp
    source_code/normalPlain.h
)

## Create the Executable  
add_executable(theExecutable
    ${CPP_H_FILES}
)

## Set Startup project
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT theExecutable)

find_package(OpenCL REQUIRED)

target_include_directories(theExecutable
    PUBLIC
        source_code
        ${OpenCL_INCLUDE_DIRS}
)


target_link_libraries(theExecutable
    OpenCL::OpenCL
)