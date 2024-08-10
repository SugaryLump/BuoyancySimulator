#include <GL/glew.h>
#include <GL/freeglut.h>

#include <vector>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

using namespace std;

class Voxels {
    private:
        vector<vector<vector<bool>>> values;
        float voxelLength;
    public:
        Voxels() = default;
        Voxels(float voxelLength, vector<vec3> boundingBox, tinyobj::attrib_t attrib, tinyobj::shape_t shape);
        void Voxels::VoxelizeVertices(vector<vec3> vertices, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner);
        void Voxels::VoxelizeEdges(vector<vector<vec3>> edges, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner);
        void Voxels::VoxelizePlanes(vec3 planeVertex, vec3 normal, float t, vector<vec3> voxelIndexBounds, vec3 minCorner);
        vector<vector<vector<bool>>> GetValues();
        float GetVoxelLength();
};