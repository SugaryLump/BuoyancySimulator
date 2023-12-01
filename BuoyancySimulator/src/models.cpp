#include <GL/glew.h>
#include <GL/freeglut.h>
#include <tiny_obj_loader.h>
#include <string>
#include <iostream>

using namespace std;

void generateVAOFromOBJ(GLuint vao, string objFilename) {
    // Load .obj
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    string err;

    bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, objFilename.c_str());

    if (!err.empty()) {
    cerr << err << endl;
    }

    if (!ret) {
    exit(1);
    }


    // Generate buffers and buffer data
    GLuint vbo, ibo;


    GLushort indices[shapes[0].mesh.indices.size()];
    for (int i = 0; i < shapes[0].mesh.indices.size(); i++) {
        //indices[i] = (short)shapes[0].mesh.indices[i];
    }
    glGenBuffers(1, &ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,)
}