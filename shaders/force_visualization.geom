#version 450

layout(triangles) in;
layout(line_strip, max_vertices=4) out;

layout(std430, binding = 1) buffer boatPositionsSSBO{
    vec4 boatPositions[];
};

layout (std430, binding = 2) buffer forcesSSBO {
    vec4 forces[];
};

layout(std430, binding = 6) buffer boatAngularPositionsSSBO {
    vec4 boatAngularPositions[];
};

layout (std430, binding = 8) buffer torquesSSBO {
    vec4 torques[];
};

uniform mat4 m_scale_rotation;
uniform mat4 m_translation;
uniform mat4 m_vp;

uniform uint boatIndex;

out vec4 colorIn;

// Rotate a vector using a given angular position
vec3 rotate(vec3 vector, vec3 angularPosition) {
    float angle = length(angularPosition);
    float cosAngle = cos(angle);
    float sinAngle = sin(angle);
    
    vec3 e = (angle == 0) ? vec3(0) : angularPosition / angle;
    return vector + (sinAngle) * (cross(e, vector)) + (1 - cosAngle) * (cross(e, cross(e, vector))); 
}

// Calculates the centroid coordinates of a triangle composed of these vertices
vec3 triangleCentroid(vec3 A, vec3 B, vec3 C) {
    return (A + B + C) / 3.0;
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

void main() {
    vec3 A = applyBoatTransforms(gl_in[0].gl_Position);
    vec3 B = applyBoatTransforms(gl_in[1].gl_Position);
    vec3 C = applyBoatTransforms(gl_in[2].gl_Position);

    vec3 center = triangleCentroid(A, B, C);
    vec3 totalTriangleForce = forces[gl_PrimitiveIDIn * 3].xyz + forces[gl_PrimitiveIDIn * 3 + 1].xyz + forces[gl_PrimitiveIDIn * 3 + 2].xyz;

    vec3 F = center + totalTriangleForce / 100000.0;
    if (totalTriangleForce.y > 0) {
        colorIn = vec4(0.0, 1.0, 0.0, 1.0);
    }
    else {
        colorIn = vec4(1.0, 0.0, 0.0, 1.0);
    }

    vec4 V1 = m_vp * vec4(center.xyz, 1.0);
    gl_Position = V1;
    EmitVertex();

    // vec4 V2 = m_vp * (vec4(normal(A, B, C), 1.0) + vec4(center.xyz, 0.0));
    vec4 V2 = m_vp * vec4(F, 1.0);
    gl_Position = V2;
    EmitVertex();

    EndPrimitive();
}