#include "models.hpp"

#include <string>
#include <vec3.hpp>
#include <mat3x3.hpp>

#pragma once
#define DEFAULT_MASS 400
#define DEFAULT_CENTER_OF_MASS vec3(0.0)
#define DEFAULT_INERTIA_MODIFIER mat3(0.003614264953);

using namespace glm;

class Buoyant {
    private:
        Model model;
        float mass;
        vec3 centerOfMass;
        mat3 inertiaModifier;
        vec3 propulsionPointOfApplication;
        vec3 propulsionForce;

    public:
        Buoyant(string boatFileName = "models/boat.boat", bool debug = false);
        Model GetModel();
        float GetMass();
        mat3 GetInertiaModifier();
        vec3 GetCenterOfMass();
        vec3 GetPropulsionPointOfApplication();
        vec3 GetPropulsionForce();
};