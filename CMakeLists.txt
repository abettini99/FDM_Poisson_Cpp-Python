cmake_minimum_required(VERSION 3.22.1)    # CMake version check
set(CMAKE_CXX_STANDARD 23)                # Enable C++23 standard
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN hidden)

## ============== ##
## Define Project ##
## ============== ##
set(PROJECT "PoissonExample")
project(PROJECT)

## ================= ##
## Create Executable ##
## ================= ##
add_executable(${PROJECT} main.cpp)

target_sources(${PROJECT}
    # no need to add headers here, only sources are required
    PRIVATE
        ${PROJECT_SOURCE_DIR}/src/insertCoefficients.cpp
)

target_include_directories(${PROJECT} 
    PRIVATE
        # where the project itself will look for internal headers
        ${PROJECT_SOURCE_DIR}/src/
    PUBLIC
        # where the project will look for public headers
        ${PROJECT_SOURCE_DIR}/external/eigen/
)

## ================= ##
## Rerout Executable ##
## ================= ##
set_target_properties(${PROJECT}
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
