#version 450

layout(triangles) in;
layout(triangle_strip, max_vertices=9) out;

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

uniform int boatIndex;
uniform float boatMass;
uniform float boatTotalArea;
uniform float boatLength;
uniform vec3 boatCenterOfMass;

uniform float waveHeight;
uniform sampler2D oceanHeightmap;

uniform mat4 m_vp;
uniform mat4 m_scale_rotation;
uniform mat4 m_translation;

out vec4 colorIn;

#define G 9.8
#define WATER_DENSITY 997.0
#define VISCOSITY_VARIATION 0.000001

// Returns a texture coordinate by normalizing a 2D vector
// within x and z [-32, 32]
vec2 heightmapCoordinate(float x, float z) {
    return vec2((x + 32)/64, (z + 32)/64);
}

float waveHeightAtPoint(float x, float z) {
    return texture(oceanHeightmap, heightmapCoordinate(x, z)).r * waveHeight;
}

// Returns a vector that is normal to 3 points
vec3 normal(vec3 A, vec3 B, vec3 C) {
    vec3 U = B - A;
    vec3 V = C - A;
    vec3 N = cross(U, V);

    if (N != vec3(0.0)) {
        N = normalize(N);
    }
    return N;
}

// Returns a plane defined by 3 points
//vec4 plane(vec3 A, vec3 B, vec3 C) {
//    vec3 normal = normal(A, B, C);
//    float delta = -dot(normal, A);
//    return vec4(normal, delta);
//}

//Returns a plane vector that roughly represents an ocean tangent plane
//for a given vertex's "horizontal" coordinates
//(x and z, since y representes altitude in this simulation)
//vec4 oceanPlane(float x, float z) {
//    float xmin = floor(x);
//    float xmax = xmin + 1;
//    float zmin = floor(z);
//    float zmax = zmin + 1;
//
//    vec3 planeA, planeB, planeC;
//
//    if (fract(x) > fract(z)) {
//        planeA = vec3(xmin, waveHeightAtPoint(xmin, zmin), zmin);
//        planeB = vec3(xmin, waveHeightAtPoint(xmin, zmax), zmax);
//        planeC = vec3(xmax, waveHeightAtPoint(xmax, zmin), zmin);
//    }
//    else {
//        planeA = vec3(xmax, waveHeightAtPoint(xmax, zmax), zmax);
//        planeB = vec3(xmax, waveHeightAtPoint(xmax, zmin), zmin);
//        planeC = vec3(xmax, waveHeightAtPoint(xmax, zmax), zmax);
//    }
//
//    return plane(planeA, planeB, planeC);
//}

// Calculate a point's signed distance to the ocean surface (negative = below surface)
// !!!Really, I have no idea why this original commented method is half as complex as it is.
// Newer version seems to work just fine and is infinitely simpler.!!!
//float surfaceDistance(vec3 vertex) {
//    vec4 pl = oceanPlane(vertex.x, vertex.z);
//
//    float oceanAltitude = -((vertex.x * pl.x) + (vertex.z * pl.z) + pl.w) / pl.y;
//    return vertex.y - oceanAltitude;
//}

// Calculate a point's signed distance to the ocean surface (negative = below surface)
float surfaceDistance(vec3 vertex) {
    return vertex.y - waveHeightAtPoint(vertex.x, vertex.z);
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

    mat3 rotationMatrix = mat3(xx + (1 - xx) * cosAngle,    xy * (1 - cosAngle) + z * sinAngle,   xz * (1 - cosAngle) - y * sinAngle,
                               xy * (1 - cosAngle) - z * sinAngle, yy + (1 - yy) * cosAngle,      yz * (1 - cosAngle) + x * sinAngle,
                               xz * (1 - cosAngle) + y * sinAngle, yz * (1 - cosAngle) - x * sinAngle,   zz + (1 - zz) * cosAngle);
    
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

vec3 applyBoatTransforms(vec4 vertex, mat4 m_model) {
    return vec3(m_model * vertex);
}

// Transforms an #workingTriangles amount of vertices (must be multiple of 3)
// in the vertex array transforming them to projection space and then emits
// the appropriate amount of triangles 
void transformAndEmitVertices(vec3[9] vertices, int workingTriangles) {
    for (int i = 0; i < workingTriangles; i++) {
        for (int v = 0; v < 3; v++) {
            vec4 vertex = m_vp * vec4(vertices[i * 3 + v], 1.0);
            //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
            //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}

// Triangle orientation check
bool compHand(vec3 A, vec3 B, vec3 C, vec3 H, vec3 M, vec3 L) {
    vec3 V1 = cross(B - A, C - A);
    vec3 V2 = cross(M - H, L - H);

    return (V1.x > 0) == (V2.x > 0) && (V1.y > 0) == (V2.y > 0) && (V1.z > 0) == (V2.z > 0);
}

// Calculates the centroid coordinates of a triangle composed of these vertices
vec3 triangleCentroid(vec3 A, vec3 B, vec3 C) {
    return (A + B + C) / 3.0;
}

// Calculates the area of a triangle composed of these vertices
float triangleArea(vec3 A, vec3 B, vec3 C) {
    return length(cross(B - A, C - A)) / 2;
}

// Calculates the buoyancy force exerted on a triangle composed of these vertices
vec3 buoyancyForce(vec3 tCentroid, vec3 tNormal, float tArea) {
    float hCenter = surfaceDistance(tCentroid);
    float force_y = (G * hCenter * WATER_DENSITY * tArea * tNormal).y;

    return vec3(0.0, force_y, 0.0);
}

// Calculates this boat's resistance coefficient
float resistanceCoefficient() {
    float reynoldsNumber = length(boatVelocities[boatIndex]) * boatLength / VISCOSITY_VARIATION;

    return 0.075 / pow(((log(reynoldsNumber) / log(10)) - 2.0), 2.0);
}

// Calculates a triangle's velocity
vec3 triangleVelocity(vec3 tCentroid) {
    vec3 vB = boatVelocities[boatIndex].xyz;
    vec3 oB = boatAngularVelocities[boatIndex].xyz;
    vec3 rBA = tCentroid - boatPositions[boatIndex].xyz;
    
    return vB + cross(oB, rBA);
}

// Calculates this triangle's viscous water resistance
// !!!The original function for this is frankly baffling. Uncountable
// expensive recalculations, completely superfluous operations that mean
// and achieve absolutely nothing. Literally calculating  the magnitude of
// a normal vector (which we KNOW is always 1 or 0) and then setting it to 1 anyway
// in case the normal is 0...
// Also the original function doesn't seem to respect its own formula,
// where it's specified that vF needs to be squared and not multiplied by it's
// length. The paper for this has it like the implementation, so I'm not sure what
// changed here.
// This is pretty modified, but if stuff is wonky, RECHECK HERE 1ST!!!
vec3 viscousWaterResistance(vec3 tNormal, float tArea, vec3 tVel, float tVelMagnitude, float resC) {
    vec3 vTangent = cross(tNormal, cross(tVel, tNormal));
    float vTangentMagnitude = length(vTangent);
    float mVel = vTangentMagnitude == 0.0 ? 1.0 : 1.0 / vTangentMagnitude;
    vec3 tangentialDir = mVel * vTangent;
    vec3 vF = tVelMagnitude * tangentialDir;
    float vFMagnitude = length(vF);
    
    return 0.5 * WATER_DENSITY * resC * tArea * vFMagnitude * vF;
}

// Calculates this triangle's cos-theta (unsure what this is)
float triangleCosTheta(vec3 tVel, vec3 tNormal) {
    vec3 normalizedTriangleVelocity = tVel == vec3(0.0) ? tVel : normalize(tVel);
    return dot(normalizedTriangleVelocity, tNormal);
}

// Calculates this triangle's pressure drag force
// !!!Again... what was going on here?? More expensive recalculations, more
// superfluous code that achieves nothing (in this case, literally,
// *programatically* nothing, even),
// checking if a power of two is negative...
// Another suspect to check if stuff is wonky!!!
vec3 pressureDragForce(vec3 tNormal, float tArea, float tVelMagnitude, float tCosTheta) {
    float C_D1, C_D2, f_P;

    // This is like this for rudimentar parametrization
    if (tCosTheta > 0) {
        C_D1 = 1;
        C_D2 = 1;
        f_P = 0.5;
    }
    else {
        C_D1 = 1;
        C_D2 = 1;
        f_P = 0.5;
    }

    return -(C_D1 * tVelMagnitude +  C_D2 * (tVelMagnitude * tVelMagnitude)) * tArea * pow(tCosTheta, f_P) * tNormal;
}

// Calculates the slamming force on this triangle
vec3 slammingForce(float tArea, float tCosTheta) {
    float boatArea = boatTotalArea;

    if (boatArea == 0) {
        return vec3(0.0);
    }

    return ((2 * boatMass * (tCosTheta < 0 ? 0 : tArea * tCosTheta)) / boatArea) * boatVelocities[boatIndex].xyz;
}

void setTriangleForceAndTorque(inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3[9] vertices, int index, float resC, vec3 worldCenterOfMass) {
    vec3 tNormal = normal(vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]);
    vec3 tCentroid = triangleCentroid(vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]);
    vec3 tVelocity = triangleVelocity(tCentroid);
    float tVelocityMagnitude = length(tVelocity);
    float tArea = triangleArea(vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]);
    float tCosTheta = triangleCosTheta(tVelocity, tNormal);

    forcesArray[index] = buoyancyForce(tCentroid, tNormal, tArea);
                         //viscousWaterResistance(tNormal, tArea, tVelocity, tVelocityMagnitude, resC);
                         //pressureDragForce(tNormal, tArea, tVelocityMagnitude, tCosTheta) +
                         //slammingForce(tArea, tCosTheta);

    // Using boat position as a center of mass is EXTREMELY sus
    // Give this a closer look eventually...
    torquesArray[index] = cross(tCentroid - worldCenterOfMass, forcesArray[index]);
}

void main() {
    vec3 A = applyBoatTransforms(gl_in[0].gl_Position);
    vec3 B = applyBoatTransforms(gl_in[1].gl_Position);
    vec3 C = applyBoatTransforms(gl_in[2].gl_Position);

    float ASurfaceDist = surfaceDistance(A);
    float BSurfaceDist = surfaceDistance(B);
    float CSurfaceDist = surfaceDistance(C);

    bool AIsSubmerged = ASurfaceDist < 0;
    bool BIsSubmerged = BSurfaceDist < 0;
    bool CIsSubmerged = CSurfaceDist < 0;


    vec3[9] vertices;


    // Auxiliary variables to cut down dozens of lines of code
    vec3[3] baseVerticesAux = {A, B, C};
    float[3] baseDistsAux = {ASurfaceDist, BSurfaceDist, CSurfaceDist};
    // Triangle generation based on how many vertices were submerged
    int totalSubmergedVertices = (AIsSubmerged ? 1 : 0) + (BIsSubmerged ? 1 : 0) + (CIsSubmerged ? 1 : 0);
    int workingTriangles, submergedTriangles;
    switch (totalSubmergedVertices) {
        case 0: {
            // No submerged triangles
            workingTriangles = 1;
            submergedTriangles = 0;

            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;
        }
            break;
        case 1: {
            // Last triangle is submerged
            workingTriangles = 3;
            submergedTriangles = 1;

            int LIndexAux = AIsSubmerged ? 0 : (BIsSubmerged ? 1 : 2);
            int MIndexAux = AIsSubmerged ? (B.y < C.y ? 1 : 2) : (BIsSubmerged ? (A.y < C.y ? 0 : 2) : (A.y < B.y ? 0 : 1));
            int HIndexAux = AIsSubmerged ? (B.y >= C.y ? 1 : 2) : (BIsSubmerged ? (A.y >= C.y ? 0 : 2) : (A.y >= B.y ? 0 : 1));

            vec3 L = baseVerticesAux[LIndexAux];
            vec3 M = baseVerticesAux[MIndexAux];
            vec3 H = baseVerticesAux[HIndexAux];
            float hL = baseDistsAux[LIndexAux];
            float hM = baseDistsAux[MIndexAux];
            float hH = baseDistsAux[HIndexAux];

            float tH = -hL / (hH - hL);
            vec3 JH = tH * (H - L) + L;
            
            float tM = -hL / (hM - hL);
            vec3 JM = tM * (M - L) + L;

            if (compHand(A, B, C, H, M, L)) {
                vertices[0] = H;
                vertices[1] = M;
                vertices[2] = JM;

                vertices[3] = JM;
                vertices[4] = JH;
                vertices[5] = H;

                vertices[6] = JM;
                vertices[7] = L;
                vertices[8] = JH;
            }
            else {
                vertices[0] = M;
                vertices[1] = H;
                vertices[2] = JM;

                vertices[3] = JH;
                vertices[4] = JM;
                vertices[5] = H;

                vertices[6] = L;
                vertices[7] = JM;
                vertices[8] = JH;
            }
        }
            break;
        case 2: {
            // Last 2 triangles are submerged (3 triangles)
            workingTriangles = 3;
            submergedTriangles = 2;
            // Vertex naming change from original implementation:
            // H -> L
            // L -> M
            // M -> H
            int LIndexAux = AIsSubmerged ? (BIsSubmerged ? 2 : 1) : 0;
            int MIndexAux = AIsSubmerged ? (BIsSubmerged ? (A.y < B.y ? 0 : 1) : (A.y < C.y ? 0 : 2)) : (B.y < C.y ? 1 : 2);
            int HIndexAux = AIsSubmerged ? (BIsSubmerged ? (A.y >= B.y ? 0 : 1) : (A.y >= C.y ? 0 : 2)) : (B.y >= C.y ? 1 : 2);

            vec3 L = baseVerticesAux[LIndexAux];
            vec3 M = baseVerticesAux[MIndexAux];
            vec3 H = baseVerticesAux[HIndexAux];
            float hL = baseDistsAux[LIndexAux];
            float hM = baseDistsAux[MIndexAux];
            float hH = baseDistsAux[HIndexAux];

            float tM = -hM / (hL - hM);
            vec3 IM = tM * (L - M) + M;

            float tH = -hH / (hL - hH);
            vec3 IH = tH * (L - H) + H;

            if (compHand(A, B, C, L, H, M)) {
                vertices[0] = L;
                vertices[1] = IH;
                vertices[2] = IM;

                vertices[3] = IH;
                vertices[4] = H;
                vertices[5] = M;

                vertices[6] = IH;
                vertices[7] = M;
                vertices[8] = IM;
            }
            else {
                vertices[0] = IH;
                vertices[1] = L;
                vertices[2] = IM;

                vertices[3] = H;
                vertices[4] = IH;
                vertices[5] = M;

                vertices[6] = M;
                vertices[7] = IH;
                vertices[8] = IM;
            }
        }
            break;
        case 3: {
            // Triangle is submerged (1 triangle)
            workingTriangles = 1;
            submergedTriangles = 1;

            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;
        }
            break;
        default: break;
    }

    vec3[3] forcesAux;
    vec3[3] torquesAux;
    forcesAux[0] = forcesAux[1] = forcesAux[2] = vec3(0.0);
    torquesAux[0] = torquesAux[1] = torquesAux[2] = vec3(0.0);
    float resC = resistanceCoefficient();
    vec3 worldCenterOfMass = applyBoatTransforms(vec4(boatCenterOfMass, 1.0));

    // case 3, where all points of the original triangle are submerged
    if (workingTriangles == 1 && submergedTriangles == 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, 0, resC, worldCenterOfMass);
    }
    // case 1 and 2, where the third generated triangle is always submerged
    else if (submergedTriangles > 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, 2, resC, worldCenterOfMass);

        // case 2, where the second generated triangle is always submerged
        if (submergedTriangles == 3) {
            setTriangleForceAndTorque(forcesAux, torquesAux, vertices, 1, resC, worldCenterOfMass);
        }
    }

    forces[gl_PrimitiveIDIn * 3] = vec4(forcesAux[0], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 1] = vec4(forcesAux[1], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 2] = vec4(forcesAux[2], 0.0);

    torques[gl_PrimitiveIDIn * 3] = vec4(torquesAux[0], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 1] = vec4(torquesAux[1], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 2] = vec4(torquesAux[2], 0.0);

    // transformAndEmitVertices(vertices, workingTriangles);
}