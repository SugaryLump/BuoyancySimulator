#include "buoyant.hpp"

#include "models.hpp"

#include <string>
#include <vec3.hpp>
#include <mat3x3.hpp>

using namespace glm;

Buoyant::Buoyant(string modelFileName,
                 float mass,
                 vec3 centerOfMass) {
    this->model = Model(modelFileName);

    this->mass = mass;
    this->centerOfMass = centerOfMass;
    this->invInertia = DEFAULT_INV_INERTIA(model.GetVolume());
}

Model Buoyant::GetModel() {
    return this->model;
}

float Buoyant::GetMass() {
    return this->mass;
}

mat3 Buoyant::GetInvInertia() {
    return this->invInertia;
}