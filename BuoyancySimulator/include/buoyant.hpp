#include "models.hpp"

#include <string>
#include <vec3.hpp>
#include <mat3x3.hpp>

#pragma once
#define DEFAULT_MASS 400
#define DEFAULT_CENTER_OF_MASS vec3(0.0)
#define DEFAULT_INERTIA_MODIFIER 0.03614264953

using namespace glm;

class Buoyant {
    private:
        Model model;
        float mass;
        vec3 centerOfMass;
        float inertiaModifier;

    public:
        Buoyant(string boatFileName = "models/boat.boat");
        Model GetModel();
        float GetMass();
        float GetInertiaModifier();
        vec3 GetCenterOfMass();
};