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

vec4 axisAngleToQuaternion(vec3 axisAngle) {
    float angle = length(axisAngle);
    if (angle == 0) return vec4(0,0,0,1);
    float sinHalfAngle = sin(angle / 2.0);
    float cosHalfAngle = cos(angle / 2.0);
    return vec4(normalize(axisAngle) * sinHalfAngle, cosHalfAngle);
}

mat3 quaternionToRotationMatrix(vec4 quat) {
    float i = quat.x;
    float j = quat.y;
    float k = quat.z;
    float r = quat.w;
    float ii = i * i;
    float ij = i * j;
    float ik = i * k;
    float ir = i * r;
    float jj = j * j;
    float jk = j * k;
    float jr = j * r;
    float kk = k * k;
    float kr = k * r;
    float rr = r * r;

    return mat3(
        -1 + 2*ii + 2*rr, 2 * (ij + kr), 2 * (ik - jr),
        2 * (ij - kr), -1 + 2*jj + 2*rr, 2 * (jk + ir),
        2 * (ik + jr), 2 * (jk - ir), -1 + 2*kk + 2*rr
    );
}

vec4 multiplyQuaternions(vec4 q1, vec4 q2) {
    mat4 auxMatrix = mat4(
        q1.w, q1.z, -q1.y, -q1.x,
        -q1.z, q1.w, q1.x, -q1.y,
        q1.y, -q1.x, q1.w, -q1.z,
        q1.x, q1.y, q1.z, q1.w
    );
    return auxMatrix * q2;
}

// Apply model transforms to trasnform vertex to world space
vec3 applyBoatTransforms(vec4 vertex) {
    vec3 worldVertex = vec3(vertex);
    // 1. & 2. Apply model scaling and then rotation
    worldVertex = vec3(m_scale_rotation * vec4(worldVertex, 1.0));
    // 3. Apply center of mass translation
    // worldVertex = worldVertex - boatCenterOfMass;
    // 4. Apply boat rotation
    worldVertex = quaternionToRotationMatrix(boatAngularPositions[boatIndex]) * worldVertex;
    // 5. Apply model translation
    worldVertex = vec3(m_translation * vec4(worldVertex, 1.0));
    // 6. Apply boat translation
    worldVertex = worldVertex + vec3(boatPositions[boatIndex]);

    return worldVertex;
}

mat3 rotateTensor (mat3 tensor, mat3 rotationMatrix) {
    return rotationMatrix * tensor * transpose(rotationMatrix);
}

vec3 nextVelocity(float deltaTime) {
    vec3 acceleration = 1.0 / boatMass * forces[0].xyz - vec3(0, 9.81, 0);
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

vec3 nextAngularVelocity(float deltaTime, mat3 inverseInertiaTensor) {
    vec3 acceleration = inverseInertiaTensor * torques[0].xyz;
    return boatAngularVelocities[boatIndex].xyz + acceleration * deltaTime;
}

// !!!Same issue as before!!!
vec4 nextAngularPosition(vec3 angularVelocity, float deltaTime) {
    vec4 velocityQuaternion = axisAngleToQuaternion(angularVelocity * deltaTime);
    return multiplyQuaternions(velocityQuaternion, boatAngularPositions[boatIndex]);
}

void main() {
    float deltaTimeSeconds = float(deltaTime) / 1000.0;

    mat3 boatRotationMatrix = quaternionToRotationMatrix(boatAngularPositions[boatIndex]);

    mat3 boatInverseInertiaTensor = inverse(rotateTensor(boatInertiaModifier, boatRotationMatrix));

    vec3 newVelocity = nextVelocity(deltaTimeSeconds);
    vec3 newPosition = nextPosition(newVelocity, deltaTimeSeconds);
    vec3 newAngularVelocity = nextAngularVelocity(deltaTimeSeconds, boatInverseInertiaTensor);
    vec4 newAngularPosition = nextAngularPosition(newAngularVelocity, deltaTimeSeconds);
    
    boatVelocities[boatIndex] = vec4(newVelocity, 0.0);
    boatPositions[boatIndex] = vec4(newPosition, 0.0);
    boatAngularVelocities[boatIndex] = vec4(newAngularVelocity, 0.0);
    boatAngularPositions[boatIndex] = newAngularPosition;

    // gl_Position = m_vp * m_model * vec4(rotate(vertex, newAngularPosition) + newPosition, 1.0);
    gl_Position = m_vp *  vec4(applyBoatTransforms(vec4(vertex, 1.0)), 1.0);
}