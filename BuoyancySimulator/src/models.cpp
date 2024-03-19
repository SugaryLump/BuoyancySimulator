#include <GL/glew.h>
#include <GL/freeglut.h>

#include "models.hpp"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <string>
#include <iostream>
#include <unordered_map>
#include <vec3.hpp>
#include <mat4x4.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <gtx/hash.hpp>
#include <vector>
#include <memory>
#include <cmath>

using namespace std;
using namespace glm;

Model::Model(string objFilename, float volume, vec3 worldPosition)
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
    glGenVertexArrays(1, vao.get());
    glBindVertexArray(*vao);

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
    this->translation = mat4({{1, 0, 0, worldPosition.x}, {0, 0, 0, worldPosition.y}, {0, 0, 0, worldPosition.z}, {0, 0, 0, 1}});
    this->rotation = mat4(1.0f);
    this->scale = mat4(1.0f);
    this->modelMatrix = mat4(1.0f);

    cout << "Finished loading model " << objFilename << endl;

    // Calculate dimensions
    this->totalArea = 0;
    float xmin, xmax, zmin, zmax;
    xmin = xmax = vertices[0];
    zmin = zmax = vertices[2];

    for (int i = 0; i < shapes[0].mesh.indices.size(); i += 3) {
        vec3 A = vec3(vertices[indices[i] * 8], vertices[indices[i] * 8 + 1], vertices[indices[i] * 8 + 2]);
        vec3 B = vec3(vertices[indices[i + 1] * 8], vertices[indices[i + 1] * 8 + 1], vertices[indices[i + 1] * 8 + 2]);
        vec3 C = vec3(vertices[indices[i + 2] * 8], vertices[indices[i + 2] * 8 + 1], vertices[indices[i + 2] * 8 + 2]);
        this->totalArea += glm::length(cross(B - A, C - A)) / 2;

        xmin = std::min(xmin, std::min(A.x, std::min(B.x, C.x)));
        xmax = std::max(xmax, std::max(A.x, std::max(B.x, C.x)));
        zmin = std::min(zmin, std::min(A.z, std::min(B.z, C.z)));
        zmax = std::max(zmax, std::max(A.z, std::max(B.z, C.z)));
    }

    this->length = std::max(xmax - xmin, zmax - zmin);
    this->volume = volume;

    cout << "Finished calculating " << objFilename << " dimensions (A = " << this->totalArea << "; V = " << this->volume << ")" << endl;
}

mat4 Model::GetModelMatrix() {
    return this->modelMatrix;
}

void Model::SetScale(vec3 scale) {
    this->scale = mat4(scale.x,   0.0,       0.0,      0.0,
                       0.0,       scale.y,   0.0,      0.0,
                       0.0,       0.0,       scale.z,  0.0,
                       0.0,       0.0,       0.0,      1.0);
    UpdateModelMatrix();
}

void Model::UpdateModelMatrix() {
    this->modelMatrix = this->translation * this->rotation * this->scale;
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