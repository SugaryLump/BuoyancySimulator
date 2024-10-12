#version 450

layout(triangles) in;
layout(line_strip, max_vertices=18) out;

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

layout (std430, binding = 9) buffer boatOldTriangleVelocitiesSSBO {
    vec4 boatOldTriangleVelocities[];
};

layout (std430, binding = 10) buffer boatOldSubmergedAreasSSBO {
    float boatOldSubmergedAreas[];
};

uniform uint maxTriangles;
uniform uint boatIndex;
uniform float boatMass;
uniform float boatTotalArea;
uniform float boatLength;
uniform vec3 boatCenterOfMass;
uniform uint oldDeltaTime;

uniform float waveHeight;
uniform sampler2D oceanHeightmap;

uniform mat4 m_vp;
uniform mat4 m_scale_rotation;
uniform mat4 m_translation;

out vec4 colorIn;

#define G 9.8
#define WATER_DENSITY 997.0
#define WATER_VISCOSITY 1.0798 // This is in CENTISTOKES!!! christ

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

// Calculate a point's signed distance to the ocean surface (negative = below surface)
float surfaceDistance(vec3 vertex) {
    return vertex.y - waveHeightAtPoint(vertex.x, vertex.z);
}

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

// Splits the triangle into two triangles with horizontal bases by
// finding a point along the line defined by highest and lowest
// vertices and connecting it to the middle point
// Places apexes of new triangles as the first element of the arrays
// The whole sorting algorithm here can be avoided if I just create my triangles
// right from the get-go. Should be changed eventually...
void horizontalBaseSplit(vec3[3] triangle, out vec3[2][3] trianglesOut) {
    // High, Middle, Middle New, Low
    int H, M, L;
    vec3 MN;

    H = 0;
    if (triangle[1].y > triangle[H].y) {
        M = 0;
        H = 1;
    }
    else {
        M = 1;
    }

    if (triangle[2].y > triangle[H].y) {
        L = M;
        M = H;
        H = 2;
    }
    else if (triangle[2].y > triangle[M].y) {
        L = M;
        M = 2;
    }
    else {
        L = 2;
    }

    vec3 upwardsTriangle[3];
    vec3 downwardsTriangle[3];

    // k is the fraction of HL that we can add to H to obtain HM
    vec3 HL = triangle[L] - triangle[H];
    float k;
    if (HL.y == 0) {
        k = 0.5;
    }
    else {
        k = (triangle[M].y - triangle[H].y) / HL.y;
    }
    MN = triangle[H] + HL * k;

    upwardsTriangle[0] = triangle[H];
    downwardsTriangle[0] = triangle[L];
    if (M == H + 1 || M == H - 2) {
        upwardsTriangle[1] = triangle[M];
        upwardsTriangle[2] = MN;
        downwardsTriangle[1] = MN;
        downwardsTriangle[2] = triangle[M];
    }
    else {
        upwardsTriangle[1] = MN;
        upwardsTriangle[2] = triangle[M];
        downwardsTriangle[1] = triangle[M];
        downwardsTriangle[2] = MN;
    }
    trianglesOut[0] = upwardsTriangle;
    trianglesOut[1] = downwardsTriangle;
}

// Calculates the real buoyancy force point of application for a given submerged
// triangle with a horizontal base.
// ASSUMES THAT THE TRIANGLE'S APEX IS A, AND THAT ITS Y IS HIGHER THAN THE BASE'S!!
vec3 upwardTrianglePointOfApplication(vec3 A, vec3 B, vec3 C) {
    float h = A.y - C.y;
    float z0 = -surfaceDistance(A);
    float tc;
    if (z0 < 0) {
        z0 = 0;
    }
    if (h == 0) {
        tc = 0.5;
    }
    else {
        tc = (4 * z0 + 3 * h) / (6 * z0 + 4 * h);
    }

    vec3 D = (C + B) / 2;
    vec3 AD = D - A;

    return A + AD * tc;
}

// Calculates the real buoyancy force point of application for a given
// submergedtriangle with a horizontal base.
// ASSUMES THAT THE TRIANGLE'S APEX IS A, AND THAT ITS Y IS LOWER THAN THE BASE'S!!
vec3 downwardTrianglePointOfApplication(vec3 A, vec3 B, vec3 C) {
    float h = C.y - A.y;
    float z0 = -surfaceDistance(B);
    float tc;
    if (z0 < 0) {
        z0 = 0;
    }
    if (h == 0) {
        tc = 0.5;
    }
    else {
        tc = (2 * z0 + h) / (6 * z0 + 2 * h);
    }
    
    vec3 D = (C + B) / 2;
    vec3 DA = A - D;

    return D + DA * tc;
}

// Calculates the area of a triangle composed of these vertices
float triangleArea(vec3 A, vec3 B, vec3 C) {
    return length(cross(B - A, C - A)) / 2;
}

// Calculates the buoyancy force exerted on a triangle composed of these vertices
vec3 buoyancyForce(vec3 pointOfApplication, vec3 tNormal, float tArea) {
    float hCenter = surfaceDistance(pointOfApplication);
    if (hCenter >= 0) {
        return vec3(0.0);
    }
    float force_y = (G * hCenter * WATER_DENSITY * tArea * tNormal).y;

    return vec3(0.0, force_y, 0.0);
}

// Calculates this boat's resistance coefficient
float resistanceCoefficient() {
    float reynoldsNumber = length(boatVelocities[boatIndex]) * boatLength * 1000000 / WATER_VISCOSITY;
    if (reynoldsNumber == 0) {
        // Should only happen if velocity is 0
        return 0;
    }

    float divisor = pow((log(reynoldsNumber-2) / log(10)), 2.0);

    float resistanceCoefficient = 0.075 / divisor;

    // Clamping because this function tends to infinity with low speed :)
    return resistanceCoefficient > 0.2 ? 0.2 : resistanceCoefficient;
}

// Calculates a point's velocity
vec3 triangleVelocity(vec3 tCentroid, vec3 rotationPoint) {
    vec3 vB = boatVelocities[boatIndex].xyz;
    vec3 oB = boatAngularVelocities[boatIndex].xyz;
    vec3 rBA = tCentroid - rotationPoint;
    
    return vB + cross(oB, rBA);
}

// Calculates this triangle's viscous water resistance
// !!!The original function for this is strange. Uncountable
// expensive recalculations, completely superfluous operations that mean
// and achieve absolutely nothing. Literally calculating  the magnitude of
// a normal vector (which we KNOW is always 1 or 0) and then setting it to 1 anyway
// in case the normal is 0... However, Kerner's formulas for this also are incorrect
// (for example, summing a vector and a scalar...)
// This is pretty modified, but if stuff is wonky, RECHECK HERE 1ST!!!
vec3 viscousWaterResistance(vec3 tNormal, float tArea, vec3 tVel, float tVelMagnitude, float resC) {
    vec3 vTangent = tVel - ((dot(tVel, tNormal) * tNormal));
    float vTangentMagnitude = length(vTangent);
    if (vTangentMagnitude == 0) {
        // This triangle's center's velocity is parallel to the its normal
        // It is moving exactly foward, and I think it's safe to say that
        // the total resulting viscous water resistance is null, since water
        // is begin dragged equally in every direction
        // Lets just throw it away :)
        return vec3(0.0);
    }
    vec3 tangentDirection = -vTangent / vTangentMagnitude;
    vec3 vF = tVelMagnitude * tangentDirection;
    
    return 0.5 * WATER_DENSITY * resC * tArea * tVelMagnitude * vF;
}

// Calculates this triangle's cos-theta
float triangleVelocityNormalCos(vec3 tVel, vec3 tNormal) {
    vec3 normalizedTriangleVelocity = tVel == vec3(0.0) ? tVel : normalize(tVel);
    return dot(normalizedTriangleVelocity, tNormal);
}

// Calculates this triangle's pressure drag force
vec3 pressureDragForce(vec3 tNormal, float tArea, float tVelMagnitude, float tCosTheta) {
    float C_D1, C_D2, f_P, speedTerm;
    int orientation;

    // This is like this for rudimentar parametrization
    if (tCosTheta > 0) {
        C_D1 = 1000;
        C_D2 = 400;
        f_P = 0.5;
        speedTerm = tVelMagnitude / 1;
        orientation = -1;
    }
    else {
        C_D1 = 1000;
        C_D2 = 150;
        f_P = 0.5;
        speedTerm = tVelMagnitude / 1;
        orientation = 1;
        tCosTheta = -tCosTheta;
    }
    // return orientation * (C_D1 * speedTerm +  C_D2 * pow(speedTerm, 2)) * tArea * pow(tCosTheta, f_P) * tNormal;
    return orientation * tNormal * (C_D1 * speedTerm + C_D2 * pow(speedTerm, 2)) * tArea * pow(tCosTheta, f_P);
}

void transformAndEmitVertices(vec3 A, vec3 B, vec3 C) {
    colorIn = vec4(1.0, 0.0, 0.0 , 1.0);
    vec4 vertex = m_vp * vec4(A, 1.0);
    gl_Position = vertex;
    EmitVertex();
    vertex = m_vp * vec4(B, 1.0);
    gl_Position = vertex;
    EmitVertex();
    vertex = m_vp * vec4(C, 1.0);
    gl_Position = vertex;
    EmitVertex();
    EndPrimitive();
}

void transformAndEmitVertices(vec3 A, vec3 B, vec3 C, float multiplier) {
    colorIn = vec4(multiplier > 0.1 ? 1 : 0, 0, 0, 1.0);
    vec4 vertex = m_vp * vec4(A, 1.0);
    gl_Position = vertex;
    EmitVertex();
    vertex = m_vp * vec4(B, 1.0);
    gl_Position = vertex;
    EmitVertex();
    vertex = m_vp * vec4(C, 1.0);
    gl_Position = vertex;
    EmitVertex();
    EndPrimitive();
}

void transformAndEmitVertices(vec3 origin, vec3 force) {
    colorIn = vec4(vec3(1), 1.0);
    gl_Position = m_vp * vec4(origin, 1.0);
    EmitVertex();
    gl_Position = m_vp * vec4(origin + force, 1.0);
    EmitVertex();
    EndPrimitive();
}

void setTriangleForceAndTorqueOLD(inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3[9] vertices, vec3 tNormal, int index, float resC, vec3 worldCenterOfMass, vec3 originPoint) {
    vec3 force = vec3(0.0);

    // Triangle splitting
    vec3 tCentroid = triangleCentroid(vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]);
    float tArea = triangleArea(vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]);
    vec3 tVelocity = triangleVelocity(tCentroid, originPoint);
    float tVelocityMagnitude = length(tVelocity);
    float tCosVelocityNormal = triangleVelocityNormalCos(tVelocity, tNormal);

    // Buoyancy Force
    force += buoyancyForce(tCentroid, tNormal, tArea);

    // Viscous Water Resistance
    force += viscousWaterResistance(tNormal, tArea, tVelocity, tVelocityMagnitude, resC);

    // Pressure Drag Force
    force += pressureDragForce(tNormal, tArea, tVelocityMagnitude, tCosVelocityNormal);

    // Torque
    vec3 torque = cross(tCentroid - worldCenterOfMass, force);

    // Final sum
    forcesArray[index] = force;
    torquesArray[index] = torque;
}

void setTriangleForceAndTorque(inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3[9] vertices, vec3 tNormal, int index, float resC, vec3 worldCenterOfMass, vec3 originPoint) {
    vec3 upForce = vec3(0.0);
    vec3 downForce = vec3(0.0);
    vec3 upTorque = vec3(0.0);
    vec3 downTorque = vec3(0.0);

    // Triangle splitting
    vec3[2][3] horBaseTriangles;
    horizontalBaseSplit(vec3[3](vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]), horBaseTriangles);
    vec3 upPointOfApplication = upwardTrianglePointOfApplication(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    vec3 downPointOfApplication = downwardTrianglePointOfApplication(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    vec3 upCentroid = triangleCentroid(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    vec3 downCentroid = triangleCentroid(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    float upArea = triangleArea(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    float downArea = triangleArea(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    vec3 upVelocity = triangleVelocity(upCentroid, originPoint);
    vec3 downVelocity = triangleVelocity(downCentroid, originPoint);
    float upVelocityMagnitude = length(upVelocity);
    float downVelocityMagnitude = length(downVelocity);
    float upCosVelocityNormal = triangleVelocityNormalCos(upVelocity, tNormal);
    float downCosVelocityNormal = triangleVelocityNormalCos(downVelocity, tNormal);

    if (upArea != 0) {
        // Buoyancy Force
        upForce += buoyancyForce(upPointOfApplication, tNormal, upArea);

        // Viscous Water Resistance
        upForce += viscousWaterResistance(tNormal, upArea, upVelocity, upVelocityMagnitude, resC);

        // Pressure Drag Force
        upForce += pressureDragForce(tNormal, upArea, upVelocityMagnitude, upCosVelocityNormal);

        // Torque
        upTorque = cross(upPointOfApplication - worldCenterOfMass, upForce);
        transformAndEmitVertices(upPointOfApplication, vec3(0, 1, 0));
    }
    if (downArea != 0) {
        // Buoyancy Force
        downForce += buoyancyForce(downPointOfApplication, tNormal, downArea);

        // Viscous Water Resistance
        downForce += viscousWaterResistance(tNormal, downArea, downVelocity, downVelocityMagnitude, resC);

    // Pressure Drag Force
    // Pressure Drag Force
    upForce += pressureDragForce(tNormal, upArea, upVelocityMagnitude, upCosVelocityNormal);
        // Pressure Drag Force
    upForce += pressureDragForce(tNormal, upArea, upVelocityMagnitude, upCosVelocityNormal);
        downForce += pressureDragForce(tNormal, downArea, downVelocityMagnitude, downCosVelocityNormal);

        // Torque
        downTorque = cross(downPointOfApplication - worldCenterOfMass, downForce);
        transformAndEmitVertices(downPointOfApplication, vec3(0, 1, 0));
    }

    // Final sum
    forcesArray[index] = upForce + downForce;
    torquesArray[index] = upTorque + downTorque;
}

void slammingForce(vec3 A, vec3 B, vec3 C, vec3 tNormal, float submergedArea, inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3 worldCenterOfMass, vec3 originPoint) {
    vec3 tCentroid = triangleCentroid(A, B, C);
    vec3 tVelocity =  triangleVelocity(tCentroid, originPoint);
    float tVelocityMagnitude = length(tVelocity);
    float tCosVelocityNormal = (triangleVelocityNormalCos(tVelocity, tNormal));

    if (submergedArea != 0 && tVelocityMagnitude != 0 && tCosVelocityNormal > 0) {
        float tOldVelocityMagnitude = length(boatOldTriangleVelocities[boatIndex * maxTriangles + gl_PrimitiveIDIn].xyz);
        float tArea = triangleArea(A, B, C);
        vec3 Fstop = -boatMass * tVelocity * 2 * submergedArea / boatTotalArea + forcesArray[0] + forcesArray[1] + forcesArray[2];
        float tOldSubmergedArea = boatOldSubmergedAreas[boatIndex * maxTriangles + gl_PrimitiveIDIn];
        vec3 gamma = (submergedArea * tVelocity - tOldSubmergedArea * boatOldTriangleVelocities[boatIndex * maxTriangles + gl_PrimitiveIDIn].xyz) / tArea;
        float gammaMax = 0.1;
        float multiplier = length(gamma) / gammaMax;
        multiplier = multiplier < 0 ? 0 : (multiplier > 1 ? 1 : multiplier);
        multiplier = pow(multiplier, 2);

        vec3 Fslam = multiplier * tCosVelocityNormal * Fstop;

        forcesArray[0] += Fslam;
        torquesArray[0] += cross(tCentroid - worldCenterOfMass, Fslam);
        
        transformAndEmitVertices(A, B, C, multiplier);
    }

    boatOldTriangleVelocities[boatIndex * maxTriangles + gl_PrimitiveIDIn] = vec4(tVelocity, 0.0);
    boatOldSubmergedAreas[boatIndex * maxTriangles + gl_PrimitiveIDIn] = submergedArea;
}

void main() {
    vec3 A = applyBoatTransforms(gl_in[0].gl_Position);
    vec3 B = applyBoatTransforms(gl_in[1].gl_Position);
    vec3 C = applyBoatTransforms(gl_in[2].gl_Position);

    vec3 tNormal = normal(A, B, C);

    float ASurfaceDist = surfaceDistance(A);
    float BSurfaceDist = surfaceDistance(B);
    float CSurfaceDist = surfaceDistance(C);

    bool AIsSubmerged = ASurfaceDist < 0;
    bool BIsSubmerged = BSurfaceDist < 0;
    bool CIsSubmerged = CSurfaceDist < 0;


    vec3[9] vertices;


    // TODO: Optimize this whole algorithm
    vec3[3] baseVerticesAux = {A, B, C};
    float[3] baseDistsAux = {ASurfaceDist, BSurfaceDist, CSurfaceDist};
    // Triangle generation based on how many vertices were submerged
    int totalSubmergedVertices = (AIsSubmerged ? 1 : 0) + (BIsSubmerged ? 1 : 0) + (CIsSubmerged ? 1 : 0);
    int workingTriangles, submergedTriangles;
    float submergedArea = 0;
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
            submergedArea += triangleArea(vertices[6], vertices[7], vertices[8]);
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
            submergedArea += triangleArea(vertices[3], vertices[4], vertices[5]);
            submergedArea += triangleArea(vertices[6], vertices[7], vertices[8]);
        }
            break;
        case 3: {
            // Triangle is submerged (1 triangle)
            workingTriangles = 1;
            submergedTriangles = 1;

            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;

            submergedArea += triangleArea(vertices[0], vertices[1], vertices[2]);
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
    vec3 originPoint = applyBoatTransforms(vec4(0, 0, 0, 1));

    // case 3, where all points of the original triangle are submerged
    if (workingTriangles == 1 && submergedTriangles == 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 0, resC, worldCenterOfMass, originPoint);
    }
    // case 1 and 2, where the third generated triangle is always submerged
    else if (workingTriangles > 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 2, resC, worldCenterOfMass, originPoint);

        // case 2, where the second generated triangle is always submerged
        if (submergedTriangles == 2) {
            setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 1, resC, worldCenterOfMass, originPoint);
        }
    }

    // slammingForce(A, B, C, tNormal, submergedArea, forcesAux, torquesAux, worldCenterOfMass, originPoint);

    forces[gl_PrimitiveIDIn * 3 + 0] = vec4(forcesAux[0], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 1] = vec4(forcesAux[1], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 2] = vec4(forcesAux[2], 0.0);

    torques[gl_PrimitiveIDIn * 3] = vec4(torquesAux[0], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 1] = vec4(torquesAux[1], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 2] = vec4(torquesAux[2], 0.0);
}