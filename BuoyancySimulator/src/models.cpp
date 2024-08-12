#include <GL/glew.h>
#include <GL/freeglut.h>

#include "models.hpp"


#include "tiny_obj_loader.h"
#include <string>
#include <iostream>
#include <unordered_map>
#include <vec3.hpp>
#include <mat4x4.hpp>

#include <gtx/hash.hpp>
#include <vector>
#include <memory>
#include <cmath>

using namespace std;
using namespace glm;

Model::Model(string objFilename, float volume, vec3 worldPosition, bool voxelsDebug)
{
    // Load .obj
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    string err;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, objFilename.c_str());

    if (!err.empty())
    {
        cerr << err << endl;
    }

    if (!ret)
    {
        exit(1);
    }

    // Restructure model data

    // This map maps the wavefront obj file's vertex indexes to a single index,
    // for IBO compatibility.
    auto index_conversion_map = unordered_map<vec3, GLuint>();
    // This vector holds our ordered vertex data (VNT)
    vector<float> vertices;
    // This array holds the mesh's shape forming indices
    GLuint indices[shapes[0].mesh.indices.size()];
    GLuint current_index = 0;

    for (int i = 0; i < shapes[0].mesh.indices.size(); i++)
    {
        tinyobj::index_t idx = shapes[0].mesh.indices[i];
        vec3 vertex_attrib_indices = {idx.vertex_index, idx.normal_index, idx.texcoord_index};

        auto emplace_result = index_conversion_map.try_emplace(
            vertex_attrib_indices,
            current_index
        );
        
        // emplace_result.first->second gets our real vertex index
        indices[i] = emplace_result.first->second;
        
        // add our new vertex data if this is a new vertex
        if (emplace_result.second) {
            vertices.push_back(attrib.vertices[3 * idx.vertex_index]);
            vertices.push_back(attrib.vertices[3 * idx.vertex_index + 1]);
            vertices.push_back(attrib.vertices[3 * idx.vertex_index + 2]);
            vertices.push_back(attrib.normals[3 * idx.normal_index]);
            vertices.push_back(attrib.normals[3 * idx.normal_index + 1]);
            vertices.push_back(attrib.normals[3 * idx.normal_index + 2]);
            vertices.push_back(attrib.texcoords[2 * idx.texcoord_index]);
            vertices.push_back(attrib.texcoords[2 * idx.texcoord_index + 1]);
            current_index++;
        }
    }

    this->minIndex = 0;
    this->maxIndex == current_index - 1;
    this->totalIndices = shapes[0].mesh.indices.size();

    // Generate IBO
    GLuint ibo, vbo;
    
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, shapes[0].mesh.indices.size() * sizeof(GLuint), indices, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Generate VBO
    int n = vertices.size();
    float vertices_aux[n];
    for (int i = 0; i < vertices.size(); i++) {
        // cout << i << "/" << vertices.size() << endl;
        vertices_aux[i] = vertices[i];
    }
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, n * sizeof(float), vertices_aux, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Generate VAO
    this->vao = make_shared<GLuint>();
    glGenVertexArrays(1, this->vao.get());
    glBindVertexArray(*(this->vao));

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(0)));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(sizeof(float)*3)));
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(sizeof(float)*6)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    
    // Unbind buffers
    glBindVertexArray(0);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
	glDisableVertexAttribArray(3);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Set matrices
    this->translation = mat4({1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {worldPosition.x, worldPosition.y, worldPosition.z, 1});
    this->rotation = mat4(1.0f);
    this->scale = mat4(1.0f);
    UpdateMatrices();

    cout << "Finished loading model " << objFilename << endl;

    // Calculate dimensions
    this->totalArea = 0;
    float xmin, xmax, zmin, zmax, ymin, ymax;
    xmin = xmax = vertices[0];
    ymin = ymax = vertices[1];
    zmin = zmax = vertices[2];

    for (int i = 0; i < shapes[0].mesh.indices.size(); i += 3) {
        vec3 A = vec3(vertices[indices[i] * 8], vertices[indices[i] * 8 + 1], vertices[indices[i] * 8 + 2]);
        vec3 B = vec3(vertices[indices[i + 1] * 8], vertices[indices[i + 1] * 8 + 1], vertices[indices[i + 1] * 8 + 2]);
        vec3 C = vec3(vertices[indices[i + 2] * 8], vertices[indices[i + 2] * 8 + 1], vertices[indices[i + 2] * 8 + 2]);
        this->totalArea += glm::length(cross(B - A, C - A)) / 2;

        xmin = std::min(xmin, std::min(A.x, std::min(B.x, C.x)));
        xmax = std::max(xmax, std::max(A.x, std::max(B.x, C.x)));
        ymin = std::min(ymin, std::min(A.y, std::min(B.y, C.y)));
        ymax = std::max(ymax, std::max(A.y, std::max(B.y, C.y)));
        zmin = std::min(zmin, std::min(A.z, std::min(B.z, C.z)));
        zmax = std::max(zmax, std::max(A.z, std::max(B.z, C.z)));
    }

    this->minBoundingBox = vector<vec3>();
    this->minBoundingBox.push_back(vec3(xmin, ymin, zmin));
    this->minBoundingBox.push_back(vec3(xmax, ymax, zmax));
    this->length = std::max(xmax - xmin, zmax - zmin);
    this->volume = volume;

    if (voxelsDebug) {
        // TODO: Clean this up, testing Voxels
        InitializeVoxelsDebug(attrib, shapes[0]);
    }
    cout << "Finished calculating " << objFilename << " dimensions (A = " << this->totalArea << "; V = " << this->volume << "; L = " << this->length << ")" << endl;
}

mat4 Model::GetTranslationMatrix() {
    return this->translation;
}

mat4 Model::GetScaleRotationMatrix() {
    return this->scaleRotationMatrix;
}

mat4 Model::GetModelMatrix() {
    return this->modelMatrix;
}

void Model::SetScale(vec3 scale) {
    this->scale = mat4(scale.x,   0.0,       0.0,      0.0,
                       0.0,       scale.y,   0.0,      0.0,
                       0.0,       0.0,       scale.z,  0.0,
                       0.0,       0.0,       0.0,      1.0);
    UpdateMatrices();
}

void Model::UpdateMatrices() {
    this->scaleRotationMatrix = this->rotation * this->scale;
    this->modelMatrix = this->translation * this->scaleRotationMatrix;
}

int Model::GetTriangleCount() {
    return this->totalIndices / 3;
}

float Model::GetTotalArea() {
    return this->totalArea;
}

float Model::GetVolume() {
    return this->volume;
}

float Model::GetLength() {
    return this->length;
}

void Model::draw() {
    glBindVertexArray(*(this->vao));
    glDrawRangeElements(GL_TRIANGLES, this->minIndex, this->maxIndex, this->totalIndices, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}

void Model::InitializeVoxelsDebug(tinyobj::attrib_t attrib, tinyobj::shape_t shape) {
    this->voxels = Voxels(0.09, this->minBoundingBox, attrib, shape);
    auto values = this->voxels.GetValues();
    auto voxelLength = this->voxels.GetVoxelLength();

    // This maps vertex position to vertex index
    auto index_conversion_map = unordered_map<vec3, GLuint>();
    // This vector holds our ordered vertex data (VNT)
    vector<vec3> vertices = vector<vec3>();
    // This array holds the mesh's shape forming indices
    vector<GLuint> indices = vector<GLuint>();
    GLuint current_index = 0;

    int xdim = values.size();
    int ydim = values[0].size();
    int zdim = values[0][0].size();

    // Create vertices
    for (int x = 0; x <= xdim; x++)
    {
        for (int y = 0; y <= ydim; y++)
        {
            for (int z = 0; z <= zdim; z++)
            {
                vertices.push_back(vec3(x * voxelLength, y * voxelLength, z * voxelLength));
            }
        }
    }

    // Create faces
    for (int x = 0; x < xdim; x++)
    {
        for (int y = 0; y < ydim; y++)
        {
            for (int z = 0; z < zdim; z++)
            {
                if (values[x][y][z]) {
                    GLuint Aindex = x * (ydim + 1) * (zdim + 1) + y * (zdim + 1) + z;
                    GLuint Bindex = Aindex + 1;
                    GLuint Cindex = Bindex + (ydim + 1) * (zdim + 1);
                    GLuint Dindex = Cindex - 1;
                    GLuint Eindex = Aindex + zdim + 1;
                    GLuint Findex = Bindex + zdim + 1;
                    GLuint Gindex = Cindex + zdim + 1;
                    GLuint Hindex = Dindex + zdim + 1;

                    // Bottom
                    indices.push_back(Aindex);
                    indices.push_back(Cindex);
                    indices.push_back(Bindex);
                    indices.push_back(Aindex);
                    indices.push_back(Dindex);
                    indices.push_back(Cindex);
                    // Top
                    indices.push_back(Eindex);
                    indices.push_back(Findex);
                    indices.push_back(Gindex);
                    indices.push_back(Eindex);
                    indices.push_back(Gindex);
                    indices.push_back(Hindex);
                    // Back
                    indices.push_back(Aindex);
                    indices.push_back(Eindex);
                    indices.push_back(Dindex);
                    indices.push_back(Eindex);
                    indices.push_back(Hindex);
                    indices.push_back(Dindex);
                    // Front
                    indices.push_back(Bindex);
                    indices.push_back(Cindex);
                    indices.push_back(Findex);
                    indices.push_back(Cindex);
                    indices.push_back(Gindex);
                    indices.push_back(Findex);
                    // Left
                    indices.push_back(Aindex);
                    indices.push_back(Bindex);
                    indices.push_back(Eindex);
                    indices.push_back(Bindex);
                    indices.push_back(Findex);
                    indices.push_back(Eindex);
                    // Right
                    indices.push_back(Dindex);
                    indices.push_back(Hindex);
                    indices.push_back(Cindex);
                    indices.push_back(Hindex);
                    indices.push_back(Gindex);
                    indices.push_back(Cindex);
                }
            }
        }
    }

    this->voxelDebugminIndex = 0;
    this->voxelDebugmaxIndex == indices.size() - 1;
    this->voxelDebugtotalIndices = indices.size();

    int n_vertices = vertices.size() * 8;
    float vertices_aux[n_vertices];
    for (int i = 0; i < vertices.size(); i++) {
        // cout << i << "/" << vertices.size() << endl;
        vertices_aux[i*8] = vertices[i].x;
        vertices_aux[i*8+1] = vertices[i].y;
        vertices_aux[i*8+2] = vertices[i].z;
        vertices_aux[i*8+3] = 0;
        vertices_aux[i*8+4] = 0;
        vertices_aux[i*8+5] = 0;
        vertices_aux[i*8+6] = 0;
        vertices_aux[i*8+7] = 0;
    }
    int n_indices = indices.size();
    GLuint indices_aux[n_indices];
    for (int i = 0; i < indices.size(); i++) {
        // cout << i << "/" << vertices.size() << endl;
        indices_aux[i] = indices[i];
    }

    // Generate IBO
    GLuint ibo, vbo;
    
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_indices * sizeof(GLuint), indices_aux, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

    // Generate VBO
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(float), vertices_aux, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Generate VAO
    this->voxelDebugvao = make_shared<GLuint>();
    glGenVertexArrays(1, this->voxelDebugvao.get());
    glBindVertexArray(*(this->voxelDebugvao));

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(0)));
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(sizeof(float)*3)));
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), ((void*)(sizeof(float)*6)));
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    
    // Unbind buffers
    glBindVertexArray(0);
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
	glDisableVertexAttribArray(3);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void Model::drawVoxelsDebug() {
    glBindVertexArray(*(this->voxelDebugvao));
    glDrawRangeElements(GL_TRIANGLES, this->voxelDebugminIndex, this->voxelDebugmaxIndex, this->voxelDebugtotalIndices, GL_UNSIGNED_INT, NULL);
    glBindVertexArray(0);
}
