#include "buoyant.hpp"

#include "models.hpp"

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vec3.hpp>
#include <mat3x3.hpp>

using namespace glm;
using namespace std;

Buoyant::Buoyant(string boatFileName) {
    this->mass = DEFAULT_MASS;
    this->centerOfMass = DEFAULT_CENTER_OF_MASS;
    this->invInertia = DEFAULT_INERTIA_MODIFIER;
    string modelFileName = "";
    float volume = 1.0;
    vec3 worldPosition = vec3(0.0, 0.0, 0.0);

    ifstream boatFile;
    boatFile.open("boatFileName");
    string entry;
    while (getline(boatFile, entry)) {
        istringstream entryStream = istringstream(entry);
        string key, value;
        getline(entryStream, key, '=');
        getline(entryStream, value);
        // TODO: Do if cases
    }



    this->model = Model(modelFileName, volume, worldPosition);
}

Model Buoyant::GetModel() {
    return this->model;
}

float Buoyant::GetMass() {
    return this->mass;
}

float Buoyant::GetInvInertia() {
    return this->invInertia;
}