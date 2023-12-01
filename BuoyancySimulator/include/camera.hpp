#pragma once

#include <vec4.hpp>

using namespace glm;

class Camera {
    private:
        float near;
        float far;
        float fov;
        vec4 position;
        vec4 lookAt;
        vec4 up;

    public:
        Camera();
        void UpdateMatrices(float ratio);
};