cmake_minimum_required(VERSION 3.10)

project(BuoyancySimulator)

# Set source files
set(SOURCES
    src/main.cpp
    src/camera.cpp
)

# Add executable
add_executable(BuoyancySimulator ${SOURCES})

# Add project include directories
message(${glm_SOURCE_DIR})
target_include_directories(BuoyancySimulator
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${freeglut_SOURCE_DIR}/include
    ${glm_SOURCE_DIR}/glm
)

# Specify library dependencies
target_link_libraries(BuoyancySimulator freeglut glm)