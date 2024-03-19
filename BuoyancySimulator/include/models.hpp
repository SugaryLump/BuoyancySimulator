#include <GL/glew.h>
#include <GL/freeglut.h>

#include <string>
#include <mat4x4.hpp>
#include <vec3.hpp>
#include <memory>

#pragma once

#define DEFAULT_VOLUME 16.918842

using namespace std;
using namespace glm;

class Model {
    private:
        mat4 translation;
        mat4 rotation;
        mat4 scale;
        mat4 modelMatrix;
        shared_ptr<GLuint> vao;
        int minIndex;
        int maxIndex;
        int totalIndices;

        float totalArea;
        float volume;
        float length;

    public:
        Model() = default;
        Model(string objFilename, float volume = 1, vec3 worldPosition = {0, 0, 0});
        mat4 GetModelMatrix();
        void SetScale(vec3 scale);
        void UpdateModelMatrix();
        int GetTriangleCount();
        float GetTotalArea();
        float GetVolume();
        float GetLength();
        void draw();
};