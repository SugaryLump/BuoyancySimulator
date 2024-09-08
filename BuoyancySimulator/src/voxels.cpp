#include <GL/glew.h>
#include <GL/freeglut.h>

#include "voxels.hpp"

#include <vector>
#include "tiny_obj_loader.h"
#include <vec3.hpp>
#include <mat4x4.hpp>
#include <cmath>
#include <string>
#include <iostream>
#include <random>

using namespace std;
using namespace glm;

#define sgn(x) (x > 0) - (x < 0)

// Voxelization algorithm taken from "Reading and voxelization of 3D models, 2013"
Voxels::Voxels(float voxelLength, vector<vec3> boundingBox, tinyobj::attrib_t attrib, tinyobj::shape_t shape) {
    int totalXVoxels = ceil((boundingBox[1].x - boundingBox[0].x) / voxelLength); 
    int totalYVoxels = ceil((boundingBox[1].y - boundingBox[0].y) / voxelLength);
    int totalZVoxels = ceil((boundingBox[1].z - boundingBox[0].z) / voxelLength);

    this->values = vector<vector<vector<bool>>>(totalXVoxels, vector<vector<bool>>(totalYVoxels, vector<bool>(totalZVoxels, false)));
    this->valuesVisited = vector<vector<vector<bool>>>(totalXVoxels, vector<vector<bool>>(totalYVoxels, vector<bool>(totalZVoxels, false)));
    this->normals = vector<vector<vector<vector<vec3>>>>(totalXVoxels, vector<vector<vector<vec3>>>(totalYVoxels, vector<vector<vec3>>(totalZVoxels, vector<vec3>())));
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
            (float)attrib.vertices[idxA.vertex_index * 3 ],
            (float)attrib.vertices[idxA.vertex_index * 3 + 1],
            (float)attrib.vertices[idxA.vertex_index * 3 + 2]
        );
        vec3 B = vec3(
            (float)attrib.vertices[idxB.vertex_index * 3 ],
            (float)attrib.vertices[idxB.vertex_index * 3 + 1],
            (float)attrib.vertices[idxB.vertex_index * 3 + 2]
        );
        vec3 C = vec3(
            (float)attrib.vertices[idxC.vertex_index * 3 ],
            (float)attrib.vertices[idxC.vertex_index * 3 + 1],
            (float)attrib.vertices[idxC.vertex_index * 3 + 2]
        );
        vec3 triangleNormal = normalize(vec3(
            (float)attrib.normals[idxA.normal_index * 3 ],
            (float)attrib.normals[idxA.normal_index * 3 + 1],
            (float)attrib.normals[idxA.normal_index * 3 + 2]
        ));
        vertices.push_back(A);
        vertices.push_back(B);
        vertices.push_back(C);
        vector<vector<vec3>> edges = vector<vector<vec3>>(3, vector<vec3>(2, vec3()));
        edges[0][0] = A;
        edges[0][1] = B;
        edges[1][0] = A;
        edges[1][1] = C;
        edges[2][0] = B;
        edges[2][1] = C;
        // Cosine of angle between voxel's diagonal and triangle normal
        float cosAlpha = abs(dot(KNormalized, triangleNormal));
        // Distance from triangle to voxel centers for voxel activation
        float t = sqrt(3) * (this->voxelLength / 2.0) * cosAlpha;
        // Triangle bounding box
        float xmin, xmax, ymin, ymax, zmin, zmax;
        xmin = std::min(A.x, std::min(B.x, C.x));
        zmin = std::min(A.z, std::min(B.z, C.z));
        ymin = std::min(A.y, std::min(B.y, C.y));
        xmax = std::max(A.x, std::max(B.x, C.x));
        zmax = std::max(A.z, std::max(B.z, C.z));
        ymax = std::max(A.y, std::max(B.y, C.y));
        vec3 minVoxelsIndices = vec3(
            std::max(0, std::min((int)((xmin - boundingBox[0].x) / voxelLength), totalXVoxels - 1)),
            std::max(0, std::min((int)((ymin - boundingBox[0].y) / voxelLength), totalXVoxels - 1)),
            std::max(0, std::min((int)((zmin - boundingBox[0].z) / voxelLength), totalXVoxels - 1))
        );
        vec3 maxVoxelsIndices = vec3(
            std::max(0, std::min((int)((xmax - boundingBox[0].x) / voxelLength), totalXVoxels - 1)),
            std::max(0, std::min((int)((ymax - boundingBox[0].y) / voxelLength), totalYVoxels - 1)),
            std::max(0, std::min((int)((zmax - boundingBox[0].z) / voxelLength), totalZVoxels - 1))
        );
        vector<vec3> voxelIndexBounds = vector<vec3>();
        voxelIndexBounds.push_back(minVoxelsIndices);
        voxelIndexBounds.push_back(maxVoxelsIndices);
        
        // Step 1 - Voxelizing Vertices
        // this->VoxelizeVertices(vertices,rC, voxelIndexBounds,boundingBox[0]);
        // Step 2 - Voxelizing Edges
        // this->VoxelizeEdges(edges, rC, voxelIndexBounds, boundingBox[0]);
        // Step 3 - Voxelizing Planes
        //this->VoxelizePlanes(vertices, triangleNormal,rC,voxelIndexBounds,boundingBox[0]);
    }
    srand((unsigned int)time(NULL));
    this->FillVolume(boundingBox[0], attrib, shape);
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

void Voxels::VoxelizePlanes(vector<vec3> triangleABC, vec3 normal, float Rc, vector<vec3> voxelIndexBounds, vec3 minCorner) {
    vec3 u = triangleABC[1] - triangleABC[0];
    vec3 v = triangleABC[2] - triangleABC[0];
    vec3 n = cross(u, v);

    float dComponent = -(n.x * triangleABC[0].x + n.y * triangleABC[0].y + n.z * triangleABC[0].z); 
    for (int x = voxelIndexBounds[0].x; x <= voxelIndexBounds[1].x; x++) {
        for (int y = voxelIndexBounds[0].y; y <= voxelIndexBounds[1].y; y++) {
            for (int z = voxelIndexBounds[0].z; z <= voxelIndexBounds[1].z; z++) {
                vec3 voxelCenter = minCorner + vec3(x * this->voxelLength, y * this->voxelLength, z * this->voxelLength) + vec3(this->voxelLength / 2);
                // check if point is above triangle using its projected baricentric coordinates
                // reference: W. Heidrich, Journal of Graphics, GPU, and Game Tools,Volume 10, Issue 3, 2005
                vec3 w = voxelCenter - triangleABC[0];
                float gamma = dot(cross(u, w), n) / dot(n, n);
                float beta = dot(cross(w, v), n) / dot(n, n);
                float alpha = 1 - gamma - beta;

                
                float d = abs(n.x * voxelCenter.x + n.y * voxelCenter.y + n.z * voxelCenter.z + dComponent) / length(n);

                if (d <= Rc && 0 <= alpha && alpha <= 1 && 0 <= beta && beta <= 1 && 0 <= gamma && gamma <= 1) {
                    this->values[x][y][z] = true;
                    this->normals[x][y][z].push_back(normal);
                }
            }
        }
    }
}

/*void Voxels::FillVolume() {
    auto startingValues = vector<vector<vector<bool>>>(this->values.size(), vector<vector<bool>>(this->values[0].size(), vector<bool>(this->values[0][0].size(), false)));
    for (int x = 0; x < startingValues.size(); x++) {
        for (int y = 0; y < startingValues[x].size(); y++) {
            for (int z = 0; z < startingValues[x][y].size(); z++) {
                startingValues[x][y][z] = bool(this->values[x][y][z]);
            }
        }
    }
    for (int x = 0; x < startingValues.size(); x++) {
        for (int y = 0; y < startingValues[x].size(); y++) {
            for (int z = 0; z < startingValues[x][y].size(); z++) {
                if (startingValues[x][y][z]) {
                    int xdir = sgn(this->normals[x][y][z][0].x);
                    int ydir = sgn(this->normals[x][y][z][0].y);
                    int zdir = sgn(this->normals[x][y][z][0].z);
                    for (int i = 1; i < this->normals[x][y][z].size(); i++) {
                        if (sgn(xdir) == sgn(this->normals[x][y][z][i].x)) {
                            xdir = sgn(this->normals[x][y][z][i].x);
                        }
                        else {
                            xdir = 0;
                        }
                        if (sgn(ydir) == sgn(this->normals[x][y][z][i].y)) {
                            ydir = sgn(this->normals[x][y][z][i].y);
                        }
                        else {
                            ydir = 0;
                        }
                        if (sgn(zdir) == sgn(this->normals[x][y][z][i].z)) {
                            zdir = sgn(this->normals[x][y][z][i].z);
                        }
                        else {
                            zdir = 0;
                        }
                    }
                    cout << "new\n";
                    this->ExpandVolume(x, y, z, xdir, ydir, zdir);
                }
            }
        }
    }
}

void Voxels::ExpandVolume(int x, int y, int z, int xdir, int ydir, int zdir) {
    if (x < this->values.size() && y < this->values[0].size() && z < this->values[0][0].size() && !this->valuesVisited[x][y][z]) { 
        cout << x << ", " << y << ", " << z << "\n";
        this->valuesVisited[x][y][z] = true;
        this->values[x][y][z] = true;
        for(int i = 0; i < this->normals[x][y][z].size(); i++) {
            if (xdir != sgn(this->normals[x][y][z][i].x)) {
                xdir = 0;
            }
            if (ydir != sgn(this->normals[x][y][z][i].y)) {
                ydir = 0;
            }
            if (zdir != sgn(this->normals[x][y][z][i].z)) {
                zdir = 0;
            }
        }
        if (xdir != 0) {
            this->ExpandVolume(x+xdir, y, z, xdir, ydir, zdir);
            if (ydir != 0) {
                this->ExpandVolume(x+xdir, y+ydir, z, xdir, ydir, zdir);
                if (zdir != 0) {
                    this->ExpandVolume(x+xdir, y+ydir, z+zdir, xdir, ydir, zdir);
                }
            }
            if (zdir != 0) {
                this->ExpandVolume(x+xdir, y, z+zdir, xdir, ydir, zdir);
            }
        }
        if (ydir != 0) {
            this->ExpandVolume(x, y+ydir, z, xdir, ydir, zdir);
            if (zdir != 0) {
                this->ExpandVolume(x, y+ydir, z+zdir, xdir, ydir, zdir);
            }
        }
        if (zdir != 0) {
            this->ExpandVolume(x, y, z+zdir, xdir, ydir, zdir);
        }
    }
}*/

void Voxels::FillVolume(vec3 minCorner, tinyobj::attrib_t attrib, tinyobj::shape_t shape) {
    std::random_device rd;
    std:mt19937 e2(rd());
    std::uniform_real_distribution<> phiDist(0, 2*M_PI);
    std::uniform_real_distribution<> costhetaDist(-1, 1);

    // Intersection algorithm from https://courses.cs.washington.edu/courses/csep557/10au/lectures/triangle_intersection.pdf
    for (int x = 0; x < this->values.size(); x++) {
        for (int y = 0; y < this->values[x].size(); y++) {
            for (int z = 0; z < this->values[x][y].size(); z++) {
                cout << "Voxel: " << x << ", " << y << ", " << z << endl;
                vec3 voxelCenter = minCorner + vec3(x * this->voxelLength, y * this->voxelLength, z * this->voxelLength) + vec3(this->voxelLength / 2);
                float shortestDist = -1;
                vec3 closestTriangleNormal;
                vec3 closestTriangleRay;
                for (int i = 0; i < 5; i++) {
                    float costheta = costhetaDist(e2);
                    float phi = phiDist(e2);
                    float theta = acos(costheta);
                    float dir_x = sin(theta) * cos(phi);
                    float dir_y = sin(theta) * sin(phi);
                    float dir_z = cos(theta);
                    vec3 dir = vec3(dir_x, dir_y, dir_z);
                    for (int i = 0; i + 2 < shape.mesh.indices.size(); i += 3) {
                        tinyobj::index_t idxA = shape.mesh.indices[i];
                        tinyobj::index_t idxB = shape.mesh.indices[i + 1];
                        tinyobj::index_t idxC = shape.mesh.indices[i + 2];

                        // Triangle vertices and edges
                        vec3 A = vec3(
                            (float)attrib.vertices[idxA.vertex_index * 3 ],
                            (float)attrib.vertices[idxA.vertex_index * 3 + 1],
                            (float)attrib.vertices[idxA.vertex_index * 3 + 2]
                        );
                        vec3 B = vec3(
                            (float)attrib.vertices[idxB.vertex_index * 3 ],
                            (float)attrib.vertices[idxB.vertex_index * 3 + 1],
                            (float)attrib.vertices[idxB.vertex_index * 3 + 2]
                        );
                        vec3 C = vec3(
                            (float)attrib.vertices[idxC.vertex_index * 3 ],
                            (float)attrib.vertices[idxC.vertex_index * 3 + 1],
                            (float)attrib.vertices[idxC.vertex_index * 3 + 2]
                        );
                        vec3 triangleNormal = normalize(vec3(
                            (float)attrib.normals[idxA.normal_index * 3 ],
                            (float)attrib.normals[idxA.normal_index * 3 + 1],
                            (float)attrib.normals[idxA.normal_index * 3 + 2]
                        ));
                        vec3 n = triangleNormal;
                        float dComponent = dot(n, A);
                        float t = (dComponent - dot(n, voxelCenter)) / (dot(n, dir));
                        vec3 Q = voxelCenter + t * dir;

                        if (abs(t) < shortestDist || shortestDist == -1) {
                            if (dot(cross((B - A), (Q - A)), n) >= 0 && dot(cross((C - B), (Q - B)), n) >= 0 && dot(cross((A - C), (Q - C)), n) >= 0) {
                                shortestDist = abs(t);
                                closestTriangleNormal = n;
                                closestTriangleRay = dir * sign(t);
                            }
                        }
                    }
                }
                
                if (shortestDist != -1 && dot(closestTriangleNormal, closestTriangleRay) > 0) {
                    this->values[x][y][z] = true;
                }
            }
        }
    }
}

vector<vector<vector<bool>>> Voxels::GetValues() {
    return this->values;
}

float Voxels::GetVoxelLength() {
    return this->voxelLength;
}