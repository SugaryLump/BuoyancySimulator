cmake_minimum_required(VERSION 3.10)

project(BuoyancySimulator)

# Set source files
set(SOURCES
    src/main.cpp
    src/camera.cpp
    src/models.cpp
    src/shader.cpp
    src/buoyant.cpp
)

# Add executable
add_executable(BuoyancySimulator ${SOURCES})

message(${CMAKE_BINARY_DIR}/contrib/tinyobjloader)
# Add project include directories
target_include_directories(BuoyancySimulator
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${freeglut_SOURCE_DIR}/include
    ${glm_SOURCE_DIR}/glm
    ${CMAKE_BINARY_DIR}/contrib/tinyobjloader
    ${glew_SOURCE_DIR}/../../include
)

# Specify library dependencies
target_link_libraries(BuoyancySimulator freeglut glm glew)