#pragma once

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <string>
#include <mat4x4.hpp>
#include <vec3.hpp>
#include <memory>

using namespace std;
using namespace glm;

void generateVAOFromOBJ(GLuint* vao, string objFilename);

class Model {
    private:
        mat4 translation;
        mat4 rotation;
        mat4 scale;
        mat4 modelMatrix;
        shared_ptr<GLuint> vao;

    public:
        Model(string objFilename);
        mat4 GetModelMatrix();
        void draw();
};