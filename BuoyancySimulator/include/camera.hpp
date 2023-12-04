
#include <vec4.hpp>
#include <vec3.hpp>
#include <vec2.hpp>
#include <mat4x4.hpp>

#pragma once

using namespace glm;

#define MOVEMENT_SPEED 0.1f
#define CAMERA_SENSITIVITY 0.001f

class Camera {
    private:
        float near;
        float far;
        float fov;
        vec3 position;
        float pitch;
        float yaw;
        vec3 lookVector;
        vec3 rightVector;
        vec3 upVector;

        mat4 viewMatrix;
        mat4 projectionMatrix;
        mat4 viewProjectionMatrix;

        bool keyboardState[256];

    public:
        float aspectRatio;
        int windowWidth;
        int windowHeight;

        Camera();

        void Update();

        void PressKey(unsigned char key);
        void LiftKey(unsigned char key);
        void UpdatePosition();
        void ProcessMouseMotion(float x, float y);

        void UpdateVectorsAndMatrices();
        mat4 GetViewMatrix();
        mat4 GetProjectionMatrix();
        mat4 GetViewProjectionMatrix();
};