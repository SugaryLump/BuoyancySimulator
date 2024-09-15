#include "buoyant.hpp"

#include "models.hpp"
#include "voxels.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vec3.hpp>
#include <mat3x3.hpp>

using namespace glm;
using namespace std;

Buoyant::Buoyant(string boatFileName, bool debug) {
    this->mass = DEFAULT_MASS;
    this->centerOfMass = DEFAULT_CENTER_OF_MASS;
    this->inertiaModifier = DEFAULT_INERTIA_MODIFIER;
    string modelFileName = "";
    float volume = 1.0;
    vec3 worldPosition = vec3(0.0, 0.0, 0.0);

    ifstream boatFile;
    boatFile.open(boatFileName);
    string entry;

    bool automateParams = false;

    cout << "Reading " << boatFileName << endl;
    while (getline(boatFile, entry)) {
        istringstream entryStream = istringstream(entry);
        string key, value;
        getline(entryStream, key, '=');
        getline(entryStream, value);

        if (key.compare("modelFile") == 0) {
            size_t lastSepIndex = boatFileName.find_last_of('/');
            if (lastSepIndex == string::npos) {
                modelFileName = value;
            }
            else {
                modelFileName = boatFileName.substr(0, lastSepIndex + 1) + value;
            }
            cout << "modelFile=" << value << endl;
        }
        else if (key.compare("mass") == 0) {
            this->mass = stof(value);
            cout << "mass=" << value << endl;
        }
        else if (key.compare("worldPosition") == 0) {
            istringstream valueStream(value);
            for (int i = 0; i < 3; i++) {
                string number;
                getline(valueStream, number, ',');
                worldPosition[i] = stof(number);
            }
            cout << "worldPosition=(" << worldPosition.x << "; " << worldPosition.y << "; " << worldPosition.z << ")" << endl;
        }
        else if (key.compare("centerOfMass") == 0) {
            istringstream valueStream(value);
            for (int i = 0; i < 3; i++) {
                string number;
                getline(valueStream, number, ',');
                this->centerOfMass[i] = stof(number);
            }
            cout << "centerOfMass=(" << this->centerOfMass.x << "; " << this->centerOfMass.y << "; " << this->centerOfMass.z << ")" << endl;
        }
        else if (key.compare("volume") == 0) {
            volume = stof(value);
            cout << "volume=" << volume << endl;
        }
        else if (key.compare("inertiaModifier") == 0) {
            float inertiaMod = stof(value);
            this->inertiaModifier = mat3(inertiaModifier);
            cout << "inertiaModifier" << inertiaMod << endl;
        }
        else if (key.compare("automateParams") == 0) {
            automateParams = value.compare("true") == 0;
        }
    }
    boatFile.close();
    
    this->model = Model(modelFileName, volume, worldPosition, debug);
    if (automateParams) {
        tinyobj::attrib_t attrib;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;

        string err;

        bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &err, modelFileName.c_str());
        Voxels voxels = this->model.GenerateVoxels(attrib, shapes[0]);
        vec3 averagePosition = voxels.GetAveragePosition();
        float estimatedVolume = voxels.GetVolume();
        mat3 inertiaMod = voxels.GetInvertedInertiaTensor(this->mass, averagePosition);
        cout << "Estimated center of mass: " << averagePosition.x << ", " << averagePosition.y << ", " << averagePosition.z << endl;
        cout << "Estimated volume: " << estimatedVolume << endl;
        cout << "Estimated inverted inertia tensor:" << endl;
        cout << inertiaMod[0][0] << ", " << inertiaMod[1][0] << ", " << inertiaMod[2][0] << endl;
        cout << inertiaMod[0][1] << ", " << inertiaMod[1][1] << ", " << inertiaMod[2][1] << endl;
        cout << inertiaMod[0][2] << ", " << inertiaMod[1][2] << ", " << inertiaMod[2][2] << endl;
        this->centerOfMass = averagePosition;
        this->model.SetVolume(estimatedVolume);
        this->inertiaModifier = inertiaMod;
    }
    cout << "Done reading " << boatFileName << endl;
}

Model Buoyant::GetModel() {
    return this->model;
}

float Buoyant::GetMass() {
    return this->mass;
}

mat3 Buoyant::GetInertiaModifier() {
    return this->inertiaModifier;
}

vec3 Buoyant::GetCenterOfMass() {
    return this->centerOfMass;
}