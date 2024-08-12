#include <GL/glew.h>
#include <GL/freeglut.h>

#include <vector>
#include <mat4x4.hpp>
#include "tiny_obj_loader.h"

#pragma once

using namespace std;
using namespace glm;

class Voxels {
    private:
        vector<vector<vector<bool>>> values;
        float voxelLength;
    public:
        Voxels() = default;
        Voxels(float voxelLength, vector<vec3> boundingBox, tinyobj::attrib_t attrib, tinyobj::shape_t shape);
        void VoxelizeVertices(vector<vec3> vertices, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner);
        void VoxelizeEdges(vector<vector<vec3>> edges, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner);
        void VoxelizePlanes(vector<vec3> triangleABC, vec3 normal, float t, vector<vec3> voxelIndexBounds, vec3 minCorner);
        vector<vector<vector<bool>>> GetValues();
        float GetVoxelLength();
};