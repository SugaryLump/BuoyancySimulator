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

uniform mat4 m_model;
uniform mat4 m_vp;

uniform int boatIndex;

out vec4 colorIn;

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

    mat3 rotationMatrix = mat3(xx + (1 - xx) * cosAngle,    xy * (1 - cosAngle) + z * sinAngle,   xz * (1 - cosAngle) - y * sinAngle,
                               xy * (1 - cosAngle) - z * sinAngle, yy + (1 - yy) * cosAngle,      yz * (1 - cosAngle) + x * sinAngle,
                               xz * (1 - cosAngle) + y * sinAngle, yz * (1 - cosAngle) - x * sinAngle,   zz + (1 - zz) * cosAngle);
    
    return rotationMatrix * vector; 
}

// Calculates the centroid coordinates of a triangle composed of these vertices
vec3 triangleCentroid(vec3 A, vec3 B, vec3 C) {
    return (A + B + C) / 3.0;
}

vec3 applyBoatTransforms(vec3 vertex) {
    return rotate(vertex, boatAngularPositions[boatIndex].xyz) + boatPositions[boatIndex].xyz;
}

void main() {
    vec3 A = (m_model * gl_in[0].gl_Position).xyz;
    vec3 B = (m_model * gl_in[0].gl_Position).xyz;
    vec3 C = (m_model * gl_in[0].gl_Position).xyz;

    A = applyBoatTransforms(A);
    B = applyBoatTransforms(B);
    C = applyBoatTransforms(C);

    vec3 center = triangleCentroid(A, B, C);
    vec3 totalTriangleForce = forces[gl_PrimitiveIDIn * 3].xyz + forces[gl_PrimitiveIDIn * 3 + 1].xyz + forces[gl_PrimitiveIDIn * 3 + 2].xyz;
    vec3 totalTriangleTorque = torques[gl_PrimitiveIDIn * 3].xyz + torques[gl_PrimitiveIDIn * 3 + 1].xyz + torques[gl_PrimitiveIDIn * 3 + 2].xyz;

    vec3 F = center + totalTriangleForce / 10.0;
    vec3 T = center + totalTriangleTorque / 10.0;

    vec4 V1 = m_vp * vec4(center.xyz, 1.0);
    gl_Position = V1;
    colorIn = vec4(0.0, length(totalTriangleForce) / 1.0, 0.0, 1.0);
    EmitVertex();

    vec4 V2 = m_vp * vec4(F, 1.0);
    gl_Position = V2;
    colorIn = vec4(0.0, length(totalTriangleForce) / 1.0, 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
    
    vec4 V3 = m_vp * vec4(center, 1.0);
    gl_Position = V3;
    colorIn = vec4(0.0, 0.0, length(totalTriangleTorque) / 1.0, 1.0);
    EmitVertex();

    vec4 V4 = m_vp * vec4(T, 1.0);
    gl_Position = V4;
    colorIn = vec4(0.0, 0.0, length(totalTriangleTorque) / 1.0, 1.0);
    EmitVertex();
    EndPrimitive();
}