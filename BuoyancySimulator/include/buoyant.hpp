#include "models.hpp"

#include <string>
#include <vec3.hpp>
#include <mat3x3.hpp>

#pragma once
#define DEFAULT_MASS 7
#define DEFAULT_CENTER_OF_MASS vec3(0.0)
#define DEFAULT_INV_INERTIA(v) {{v/(63.241470*DEFAULT_MASS),0,0}, {0,v/(63.267620*DEFAULT_MASS), 0},{ 0,0,v/(63.241470*DEFAULT_MASS)}}

using namespace glm;

class Buoyant {
    private:
        Model model;
        float mass;
        vec3 centerOfMass;
        mat3 invInertia;

    public:
        Buoyant(string modelFileName = "models/boat.hpp",
                float mass = DEFAULT_MASS,
                vec3 centerOfMass = DEFAULT_CENTER_OF_MASS);
        Model GetModel();
        float GetMass();
        mat3 GetInvInertia();
};