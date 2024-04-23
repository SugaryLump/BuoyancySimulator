#version 450

layout(triangles) in;
layout(triangle_strip, max_vertices=15) out;

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

uniform int boatIndex;
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
#define WATER_VISCOSITY 0.00109

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
    if (hCenter > 0) {
        return vec3(0.0);
    }
    float force_y = (G * hCenter * WATER_DENSITY * tArea * tNormal).y;

    return vec3(0.0, force_y, 0.0);
}

// Calculates this boat's resistance coefficient
float resistanceCoefficient() {
    float reynoldsNumber = length(boatVelocities[boatIndex]) * boatLength / WATER_VISCOSITY;
    if (reynoldsNumber == 0) {
        // Should only happen if velocity is 0
        return 0;
    }
    // This divisor goes against the formula specified in Kerner's implementation
    // This is because the formula supplied there behaves very weirdly - close to 0 speed, there is an
    // extremely sharp spike in the coefficient, which really doesn't make any sense. The result was
    // that when a boat was begining to settle down, it would hit that spike briefly and get launched
    // incredibly fast by an unreasonable resistance force.
    // This function behaves very simillarly to the original, but doesn't have any spikes past v=0
    float divisor = pow(((log(reynoldsNumber+1) / log(10))), 3.0)/14;
    if (divisor == 0) {
        return 0;
    }

    return 0.075 / divisor;
}

// Calculates a triangle's velocity
vec3 triangleVelocity(vec3 tCentroid, vec3 worldCenterOfMass) {
    vec3 vB = boatVelocities[boatIndex].xyz;
    vec3 oB = boatAngularVelocities[boatIndex].xyz;
    vec3 rBA = tCentroid - worldCenterOfMass;
    
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
    vec3 vTangent = cross(tNormal, cross(tVel, tNormal));
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
    float vFMagnitude = length(vF);
    
    return 0.5 * WATER_DENSITY * resC * tArea * vFMagnitude * vF;
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
        C_D1 = 3;
        C_D2 = 3;
        f_P = 0.5;
        speedTerm = tVelMagnitude / 0.1;
        orientation = -1;
    }
    else {
        C_D1 = 3;
        C_D2 = 3;
        f_P = 0.5;
        speedTerm = tVelMagnitude / 0.1;
        orientation = 1;
        tCosTheta = -tCosTheta;
    }
    // return orientation * (C_D1 * speedTerm +  C_D2 * pow(speedTerm, 2)) * tArea * pow(tCosTheta, f_P) * tNormal;
    return orientation * tNormal * (C_D1 * speedTerm + C_D2 * pow(speedTerm, 2)) * tArea * pow(tCosTheta, f_P);
}

// Calculates the slamming force on this triangle
// The original function for this is straight up just using the wrong formula
/*vec3 slammingForce(float tArea, float tCosTheta) {
    float boatArea = boatTotalArea;

    if (boatArea == 0) {
        return vec3(0.0);
    }

    return ((2 * boatMass * (tCosTheta < 0 ? 0 : tArea * tCosTheta)) / boatArea) * boatVelocities[boatIndex].xyz;
}*/

// Transforms an #workingTriangles amount of vertices (must be multiple of 3)
// in the vertex array transforming them to projection space and then emits
// the appropriate amount of triangles 
void transformAndEmitVertices(vec3[9] vertices, int workingTriangles) {
    for (int i = 0; i < workingTriangles; i++) {
        vec3 center = triangleCentroid(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
        colorIn = vec4(surfaceDistance(center), 0.0,0.0,1.0);
        //colorIn = vec4(1.0, 0.0, 0.0, 1.0);
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

void transformAndEmitVertices(vec3[2][3] vertices) {
    for (int i = 0; i < 2; i++) {
        colorIn = vec4(i, 0.35, 1-i, 1.0);
        for (int v = 0; v < 3; v++) {
            //vec4 vertex = m_vp * (vec4(vertices[i][v], 1.0) + 0.01 * vec4(normal(vertices[i][0], vertices[i][1], vertices[i][2]), 0.0) + (i*2-1) * vec4(0, 0.05, 0, 0));
            vec4 vertex = m_vp * vec4(vertices[i][v], 1.0);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}

void transformAndEmitVertices(vec3[2][3] vertices, float upArea, float downArea, float tArea) {
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            colorIn = vec4(upArea/tArea, 0.5, 0.5, 0.8);
        }
        else {
            colorIn = vec4(downArea/tArea, 0.5, 0.5, 0.8);
        }
        if (colorIn.x > 1.001) {
            colorIn = vec4(0, 1, 0 , 1);
        }
        for (int v = 0; v < 3; v++) {
            vec4 vertex = m_vp * (vec4(vertices[i][v], 1.0) + 0.01 * vec4(normal(vertices[i][0], vertices[i][1], vertices[i][2]), 0.0));
            //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
            //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}

void transformAndEmitVertices(vec3[2][3] vertices, vec3 n) {
    for (int i = 0; i < 2; i++) {
        colorIn = abs(vec4(n/2, 1.0));
        for (int v = 0; v < 3; v++) {
            vec4 vertex = m_vp * (vec4(vertices[i][v], 1.0) + 0.01 * vec4(normal(vertices[i][0], vertices[i][1], vertices[i][2]), 0.0));
            //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
            //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}

void transformAndEmitVerticesForces(vec3[2][3] vertices, vec3 upForce, vec3 downForce) {
    for (int i = 0; i < 2; i++) {
        if (i == 0) {
            colorIn = vec4(0, 0.4, 1, length(upForce)/length(downForce)/3);
        }
        else {
            colorIn = vec4(1, 0.4, 0, length(downForce)/length(upForce)/3);
        }
        for (int v = 0; v < 3; v++) {
            vec4 vertex = m_vp * (vec4(vertices[i][v], 1.0) + 0.01 * vec4(normal(vertices[i][0], vertices[i][1], vertices[i][2]), 0.0));
            //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
            //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}

void transformAndEmitVertices(vec3[2][3] vertices, vec3 upPointOfApplication, vec3 downPointOfApplication, vec3 tCentroid) {
    for (int i = 0; i < 2; i++) {
        colorIn = vec4(1, 1, 1, 0.5);
        if (i == 0) {
            vec3 upcentroid = triangleCentroid(vertices[0][0], vertices[0][1], vertices[0][2]);
            if (length(upPointOfApplication - upcentroid) > 0.2) {
                colorIn = vec4(1, 0, 0, 1);
            }
        }
        else {
            vec3 downcentroid = triangleCentroid(vertices[1][0], vertices[1][1], vertices[1][2]);
            if (length(downPointOfApplication - downcentroid) > 0.2) {
                colorIn = vec4(0, 0, 1, 1);
            }
        }
        for (int v = 0; v < 3; v++) {
            vec4 vertex = m_vp * (vec4(vertices[i][v], 1.0) + 0.01 * vec4(normal(vertices[i][0], vertices[i][1], vertices[i][2]), 0.0));
            //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
            //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
            gl_Position = vertex;
            EmitVertex();
        }
        EndPrimitive();
    }
    
}


void transformAndEmitVertices(vec3[9] vertices, int i, bool dummy) {
    float area = triangleArea(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
    colorIn = vec4(1.0, 1.0, 1.0, 1.0);
    for (int v = 0; v < 3; v++) {
        vec4 vertex = m_vp * vec4(vertices[i*3 + v], 1.0);
        //colorIn = vec4(waveHeightAtPoint(vertices[i * 3 + v].x, vertices[i * 3 + v].z) / waveHeight);
        //colorIn = vec4(surfaceDistance(vertices[i * 3 + v]) * 10);
        gl_Position = vertex;
        EmitVertex();
    }
    EndPrimitive();    
}

void transformAndEmitVertices(vec3 A, vec3 B, vec3 C, float multiplier) {
    if (multiplier > 0.2) {
        colorIn = vec4(multiplier, 1-multiplier, 1-multiplier, 1.0);
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
}

void transformAndEmitVertices(vec3 A, vec3 B, vec3 C, float subArea, bool dummy) {
    float tArea = triangleArea(A, B, C);
    colorIn = vec4(subArea/tArea, 0.0, 0.0 , 1.0);
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

void setTriangleForceAndTorque(inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3[9] vertices, vec3 tNormal, int index, float resC, vec3 worldCenterOfMass) {
    vec3 upForce = vec3(0.0);
    vec3 downForce = vec3(0.0);

    // Triangle splitting
    vec3[2][3] horBaseTriangles;
    horizontalBaseSplit(vec3[3](vertices[index * 3], vertices[index * 3 + 1], vertices[index * 3 + 2]), horBaseTriangles);
    vec3 upPointOfApplication = upwardTrianglePointOfApplication(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    vec3 downPointOfApplication = downwardTrianglePointOfApplication(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    vec3 upCentroid = triangleCentroid(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    vec3 downCentroid = triangleCentroid(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    float upArea = triangleArea(horBaseTriangles[0][0], horBaseTriangles[0][1], horBaseTriangles[0][2]);
    float downArea = triangleArea(horBaseTriangles[1][0], horBaseTriangles[1][1], horBaseTriangles[1][2]);
    vec3 upVelocity = triangleVelocity(upCentroid, worldCenterOfMass);
    vec3 downVelocity = triangleVelocity(downCentroid, worldCenterOfMass);
    float upVelocityMagnitude = length(upVelocity);
    float downVelocityMagnitude = length(downVelocity);
    float upCosVelocityNormal = triangleVelocityNormalCos(upVelocity, tNormal);
    float downCosVelocityNormal = triangleVelocityNormalCos(downVelocity, tNormal);

    // Buoyancy Force
    upForce += buoyancyForce(upPointOfApplication, tNormal, upArea);
    downForce += buoyancyForce(downPointOfApplication, tNormal, downArea);

    // Viscous Water Resistance
    upForce += viscousWaterResistance(tNormal, upArea, upVelocity, upVelocityMagnitude, resC);
    downForce += viscousWaterResistance(tNormal, downArea, downVelocity, downVelocityMagnitude, resC);

    // Pressure Drag Force
    //upForce += pressureDragForce(tNormal, upArea, upVelocityMagnitude, upCosVelocityNormal);
    //downForce += pressureDragForce(tNormal, downArea, downVelocityMagnitude, downCosVelocityNormal);

    // Torque
    vec3 upTorque = cross(upPointOfApplication - worldCenterOfMass, upForce);
    vec3 downTorque = cross(downPointOfApplication - worldCenterOfMass, downForce);

    // Final sum
    forcesArray[index] = upForce + downForce;
    torquesArray[index] = upTorque + downTorque;
}

void slammingForce(vec3 A, vec3 B, vec3 C, vec3 tNormal, float submergedArea, inout vec3[3] forcesArray, inout vec3[3] torquesArray, vec3 worldCenterOfMass) {
    vec3 tCentroid = triangleCentroid(A, B, C);
    vec3 tVelocity =  triangleVelocity(tCentroid, worldCenterOfMass);
    float tArea = triangleArea(A, B, C);
    float tCosVelocityNormal = (triangleVelocityNormalCos(tVelocity, tNormal));
    if (submergedArea == 0 || oldDeltaTime == 0 || tVelocity == vec3(0) || tCosVelocityNormal <= 0) {
        boatOldTriangleVelocities[gl_PrimitiveIDIn] = vec4(tVelocity, 0.0);
        boatOldSubmergedAreas[gl_PrimitiveIDIn] = submergedArea;
        return;
    }
    float oldDeltaTimeSeconds = float(oldDeltaTime) / 1000.0;
    vec3 Fstop = -boatMass * tVelocity * 2 * submergedArea / boatTotalArea;
    float tOldSubmergedArea = boatOldSubmergedAreas[gl_PrimitiveIDIn];
    vec3 tOldVelocity = boatOldTriangleVelocities[gl_PrimitiveIDIn].xyz;
    float gamma = (submergedArea * length(tVelocity) - tOldSubmergedArea * length(tOldVelocity)) / (tArea * oldDeltaTimeSeconds);
    float gammaMax = length(tVelocity) / oldDeltaTimeSeconds / 5;
    float multiplier = gamma / gammaMax;
    multiplier = multiplier < 0 ? 0 : (multiplier > 1 ? 1 : multiplier);
    multiplier = pow(multiplier, 2);

    transformAndEmitVertices(A, B, C, multiplier);

    vec3 Fslam = multiplier * tCosVelocityNormal * Fstop;

    forcesArray[0] += Fslam;
    torquesArray[0] += cross(tCentroid - worldCenterOfMass, Fslam);

    boatOldTriangleVelocities[gl_PrimitiveIDIn] = vec4(tVelocity, 0.0);
    boatOldSubmergedAreas[gl_PrimitiveIDIn] = submergedArea;
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

    // case 3, where all points of the original triangle are submerged
    if (workingTriangles == 1 && submergedTriangles == 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 0, resC, worldCenterOfMass);
    }
    // case 1 and 2, where the third generated triangle is always submerged
    else if (workingTriangles > 1) {
        setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 2, resC, worldCenterOfMass);

        // case 2, where the second generated triangle is always submerged
        if (submergedTriangles == 2) {
            setTriangleForceAndTorque(forcesAux, torquesAux, vertices, tNormal, 1, resC, worldCenterOfMass);
        }
    }

    slammingForce(A, B, C, tNormal, submergedArea, forcesAux, torquesAux, worldCenterOfMass);

    forces[gl_PrimitiveIDIn * 3 + 0] = vec4(forcesAux[0], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 1] = vec4(forcesAux[1], 0.0);
    forces[gl_PrimitiveIDIn * 3 + 2] = vec4(forcesAux[2], 0.0);

    torques[gl_PrimitiveIDIn * 3] = vec4(torquesAux[0], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 1] = vec4(torquesAux[1], 0.0);
    torques[gl_PrimitiveIDIn * 3 + 2] = vec4(torquesAux[2], 0.0);

    //transformAndEmitVertices(vertices, workingTriangles);
}