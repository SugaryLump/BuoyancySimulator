#version 450

// This entire shade raises some questions with the entire SSBO approach.
// Using SSBOs that will be written and read the exact same way across every
// single invocation, as well as recalculating the exact same values,
// seems unnecessary. However, the cost of passing this data to
// and from GPU/CPU might not make it worthwhile to do any other way
// Worth investigating?

in vec3 vertex;

layout(std430, binding = 1) buffer boatPositionsSSBO{
    vec4 boatPositions[];
};

layout (std430, binding = 2) buffer forcesSSBO {
    vec4 forces[];
};

layout (std430, binding = 4) buffer boatVelocitiesSSBO {
    vec4 boatVelocities[];
};

layout (std430, binding = 5) buffer boatAngularVelocitiesSSBO{
    vec4 boatAngularVelocities[];
};

layout(std430, binding = 6) buffer boatAngularPositionsSSBO {
    vec4 boatAngularPositions[];
};

layout (std430, binding = 8) buffer torquesSSBO {
    vec4 torques[];
};

uniform uint boatIndex;
uniform float boatMass;
uniform mat3 boatInertiaModifier;

uniform uint deltaTime;

uniform mat4 m_vp;
uniform mat4 m_scale_rotation;
uniform mat4 m_translation;

vec3 nextVelocity(vec3 F, float deltaTime) {
    vec3 acceleration = 1.0 / boatMass * F;
    return boatVelocities[boatIndex].xyz + deltaTime  * acceleration;
}

// !!!Original source for this function is slowing the speed down by 50% and
// also adding the new velocity to the old velocity (even though the new
// velocity is already calculated by summing acceleration to the old
// velocity...)
// I don't know!!!
vec3 nextPosition(vec3 velocity, float deltaTime) {
    return boatPositions[boatIndex].xyz + deltaTime * velocity;
}

vec3 nextAngularVelocity(float deltaTime) {
    vec3 acceleration = boatInertiaModifier * torques[0].xyz;
    return boatAngularVelocities[boatIndex].xyz + acceleration * deltaTime;
}

// !!!Same issue as before!!!
vec3 nextAngularPosition(vec3 angularVelocity, float deltaTime) {
    return boatAngularPositions[boatIndex].xyz + angularVelocity * deltaTime;
}

// Rotate a vector using a given angular position
vec3 rotate(vec3 vector, vec3 angularPosition) {
    float angle = length(angularPosition);
    float cosAngle = cos(angle);
    float sinAngle = sin(angle);
    
    float x = (angle == 0) ? 0 : angularPosition.x / angle;
    float y = (angle == 0) ? 0 : angularPosition.y / angle;
    float z = (angle == 0) ? 0 : angularPosition.z / angle;
    
    float xx = x * x;
    float xy = x * y;
    float xz = x * z;
    float yy = y * y;
    float yz = y * z;
    float zz = z * z;

    mat3 rotationMatrix = mat3(cosAngle + xx * (1 - cosAngle),     xy * (1 - cosAngle) - z * sinAngle,   xz * (1 - cosAngle) + y * sinAngle,
                               xy * (1 - cosAngle) + z * sinAngle, cosAngle + yy * (1 - cosAngle),       yz * (1 - cosAngle) - x * sinAngle,
                               xz * (1 - cosAngle) - y * sinAngle, yz * (1 - cosAngle) + x * sinAngle,   cosAngle + zz * (1 - cosAngle));
    
    return rotationMatrix * vector; 
}

// Apply model transforms to trasnform vertex to world space
vec3 applyBoatTransforms(vec4 vertex) {
    vec3 worldVertex = vec3(vertex);
    // 1. & 2. Apply model scaling and then rotation
    worldVertex = vec3(m_scale_rotation * vec4(worldVertex, 1.0));
    // 3. Apply center of mass translation
    // worldVertex = worldVertex - boatCenterOfMass;
    // 4. Apply boat rotation
    worldVertex = rotate(worldVertex, vec3(boatAngularPositions[boatIndex]));
    // 5. Apply model translation
    worldVertex = vec3(m_translation * vec4(worldVertex, 1.0));
    // 6. Apply boat translation
    worldVertex = worldVertex + vec3(boatPositions[boatIndex]);

    return worldVertex;
}

void main() {
    vec3 gravity = boatMass * vec3(0.0, -9.8, 0.0);

    vec3 F = forces[0].xyz + gravity;

    float deltaTimeSeconds = float(deltaTime) / 1000.0;

    vec3 newVelocity = nextVelocity(F, deltaTimeSeconds);
    vec3 newPosition = nextPosition(newVelocity, deltaTimeSeconds);
    vec3 newAngularVelocity = nextAngularVelocity(deltaTimeSeconds);
    vec3 newAngularPosition = nextAngularPosition(newAngularVelocity, deltaTimeSeconds);
    
    boatVelocities[boatIndex] = vec4(newVelocity, 0.0);
    boatPositions[boatIndex] = vec4(newPosition, 0.0);
    boatAngularVelocities[boatIndex] = vec4(newAngularVelocity, 0.0);
    boatAngularPositions[boatIndex] = vec4(newAngularPosition, 0.0);

    // gl_Position = m_vp * m_model * vec4(rotate(vertex, newAngularPosition) + newPosition, 1.0);
    gl_Position = m_vp *  vec4(applyBoatTransforms(vec4(vertex, 1.0)), 1.0);
}