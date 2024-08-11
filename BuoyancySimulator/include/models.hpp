#include <GL/glew.h>
#include <GL/freeglut.h>

#include <string>
#include <mat4x4.hpp>
#include <vec3.hpp>
#include <memory>
#include "voxels.hpp"

#pragma once
#define DEFAULT_VOLUME 16.918842

using namespace std;
using namespace glm;

class Model {
    private:
        mat4 translation;
        mat4 rotation;
        mat4 scale;
        mat4 scaleRotationMatrix;
        mat4 modelMatrix;
        shared_ptr<GLuint> vao;
        int minIndex;
        int maxIndex;
        int totalIndices;

        vector<vec3> minBoundingBox;
        float totalArea;
        float volume;
        float length;

        Voxels voxels;
        shared_ptr<GLuint> voxelDebugvao;
        int voxelDebugminIndex;
        int voxelDebugmaxIndex;
        int voxelDebugtotalIndices;


    public:
        Model() = default;
        Model(string objFilename, float volume = 1, vec3 worldPosition = {0.0, 0.0, 0.0}, bool voxelsDebug = false);
        mat4 GetTranslationMatrix();
        mat4 GetScaleRotationMatrix();
        mat4 GetModelMatrix();
        void SetScale(vec3 scale);
        void UpdateMatrices();
        int GetTriangleCount();
        vector<vec3> GetMinBoundingBox();
        float GetTotalArea();
        float GetVolume();
        float GetLength();
        void draw();
        void drawVoxelsDebug();
        void InitializeVoxelsDebug(tinyobj::attrib_t attrib, tinyobj::shape_t shape);
};