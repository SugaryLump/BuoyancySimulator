#include <GL/glew.h>
#include <GL/freeglut.h>

#include <string>
#include <vector>

#pragma once

using namespace std;

class Shader {
    private:
    public:
        GLuint shaderProgram;

        Shader(string shaderName);
        void BindShader();
        void SetupUniforms();
};