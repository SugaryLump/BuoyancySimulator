#include <GL/glew.h>
#include <GL/freeglut.h>

#include "voxels.hpp"

#include <vector>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <vec3.hpp>
#include <mat4x4.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <cmath>

using namespace std;
using namespace glm;

// Voxelization algorithm taken from "Reading and voxelization of 3D models, 2013"
Voxels::Voxels(float voxelLength, vector<vec3> boundingBox, tinyobj::attrib_t attrib, tinyobj::shape_t shape) {
    this->values = vector<vector<vector<bool>>>();
    this->voxelLength = voxelLength;

    // Voxel's diagonal vector, normalized 
    vec3 KNormalized = normalize(vec3(1, 1, 1));
    // Distance from vertices and edges to voxel centers for voxel activation
    float rC = (sqrt(3) / 2) * this->voxelLength;

    for (int i = 0; i + 2 < shape.mesh.indices.size(); i += 3) {
        tinyobj::index_t idxA = shape.mesh.indices[i];
        tinyobj::index_t idxB = shape.mesh.indices[i + 1];
        tinyobj::index_t idxC = shape.mesh.indices[i + 2];

        // Triangle vertices and edges
        vector<vec3> vertices = vector<vec3>();
        vec3 A = vec3(
            (float)attrib.vertices[idxA.vertex_index],
            (float)attrib.vertices[idxA.vertex_index + 1],
            (float)attrib.vertices[idxA.vertex_index + 2]
        );
        vec3 B = vec3(
            (float)attrib.vertices[idxB.vertex_index],
            (float)attrib.vertices[idxB.vertex_index + 1],
            (float)attrib.vertices[idxB.vertex_index + 2]
        );
        vec3 C = vec3(
            (float)attrib.vertices[idxC.vertex_index],
            (float)attrib.vertices[idxC.vertex_index + 1],
            (float)attrib.vertices[idxC.vertex_index + 2]
        );
        vec3 triangleNormal = normalize(vec3(
            (float)attrib.normals[idxA.normal_index],
            (float)attrib.normals[idxA.normal_index + 1],
            (float)attrib.normals[idxA.normal_index + 2]
        ));
        vertices.push_back(A);
        vertices.push_back(B);
        vertices.push_back(C);
        vector<vector<vec3>> edges = vector<vector<vec3>>();
        edges[0][0] = A;
        edges[0][1] = B;
        edges[1][0] = A;
        edges[1][1] = C;
        edges[2][0] = B;
        edges[2][1] = C;
        // Cosine of angle between voxel's diagonal and triangle normal
        float cosAlpha = dot(KNormalized, triangleNormal);
        // Distance from triangle to voxel centers for voxel activation
        float t = sqrt(3) * (this->voxelLength / 2) * cosAlpha;
        // Triangle bounding box
        float xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = std::min(A.x, std::min(B.x, C.x));
        zmin = std::min(A.z, std::min(B.z, C.z));
        ymin = std::min(A.y, std::min(B.y, C.y));
        xmax = std::max(A.x, std::max(B.x, C.x));
        zmax = std::max(A.z, std::max(B.z, C.z));
        ymax = std::max(A.y, std::max(B.y, C.y));
        vec3 minVoxelsIndices = vec3(
            (int)((xmin - boundingBox[0].x) / voxelLength),
            (int)((ymin - boundingBox[0].y) / voxelLength),
            (int)((zmin - boundingBox[0].z) / voxelLength)
        );
        vec3 maxVoxelsIndices = vec3(
            (int)((xmin - boundingBox[0].x) / voxelLength) + 1,
            (int)((ymin - boundingBox[0].y) / voxelLength) + 1,
            (int)((zmin - boundingBox[0].z) / voxelLength) + 1
        );
        vector<vec3> voxelIndexBounds = vector<vec3>();
        voxelIndexBounds.push_back(minVoxelsIndices);
        voxelIndexBounds.push_back(maxVoxelsIndices);
        
        // Step 1 - Voxelizing Vertices
        this->VoxelizeVertices(vertices,rC, voxelIndexBounds,boundingBox[0]);
        // Step 2 - Voxelizing Edges
        this->VoxelizeEdges(edges, rC, voxelIndexBounds, boundingBox[0]);
        // Step 3 - Voxelizing Planes
        this->VoxelizePlanes(A, triangleNormal,t,voxelIndexBounds,boundingBox[0]);
    }
}

void Voxels::VoxelizeVertices(vector<vec3> vertices, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner) {
    for (int i = 0; i < vertices.size(); i++) {
        for (int x = voxelIndexBounds[0].x; x <= voxelIndexBounds[1].x; x++) {
            for (int y = voxelIndexBounds[0].y; y <= voxelIndexBounds[1].y; y++) {
                for (int z = voxelIndexBounds[0].z; z <= voxelIndexBounds[1].z; z++) {
                    vec3 voxelCenter = minCorner + vec3(x * this->voxelLength, y * this->voxelLength, z * this->voxelLength) + vec3(this->voxelLength / 2);
                    if (distance(voxelCenter, vertices[i]) <= Rc) {
                        this->values[x][y][z] = true;
                    }
                }
            }
        }
    }
}

void Voxels::VoxelizeEdges(vector<vector<vec3>> edges, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner) {
    for (int i = 0; i < edges.size(); i++) {
        for (int x = voxelIndexBounds[0].x; x <= voxelIndexBounds[1].x; x++) {
            for (int y = voxelIndexBounds[0].y; y <= voxelIndexBounds[1].y; y++) {
                for (int z = voxelIndexBounds[0].z; z <= voxelIndexBounds[1].z; z++) {
                    vec3 voxelCenter = minCorner + vec3(x * this->voxelLength, y * this->voxelLength, z * this->voxelLength) + vec3(this->voxelLength / 2);
                    vec3 AVc = voxelCenter - edges[i][0];
                    vec3 AB = edges[i][1] - edges[i][0];
                    float projectedLength = dot(AVc, AB) / length(AB);
                    float d = distance(edges[i][0] + (normalize(AB) * projectedLength), voxelCenter);
                    if (d <= Rc) {
                        this->values[x][y][z] = true;
                    }
                }
            }
        }
    }
}

void Voxels::VoxelizePlanes(vec3 planeVertex, vec3 normal, float t, vector<vec3> voxelIndexBounds, vec3 minCorner) {
    for (int x = voxelIndexBounds[0].x; x <= voxelIndexBounds[1].x; x++) {
        for (int y = voxelIndexBounds[0].y; y <= voxelIndexBounds[1].y; y++) {
            for (int z = voxelIndexBounds[0].z; z <= voxelIndexBounds[1].z; z++) {
                vec3 voxelCenter = minCorner + vec3(x * this->voxelLength, y * this->voxelLength, z * this->voxelLength) + vec3(this->voxelLength / 2);
                float dComponent = -(normal.x * planeVertex.x + normal.y * planeVertex.y + normal.z * planeVertex.z); 
                float d = length(normal.x * voxelCenter.x + normal.y * voxelCenter.y + normal.z * voxelCenter.z + dComponent) / length(normal);
                if (d <= t) {
                    this->values[x][y][z] = true;
                }
            }
        }
    }
}
