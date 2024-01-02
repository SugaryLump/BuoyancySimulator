#version 450

layout(triangles) in;
layout(triangle_strip, max_vertices=9) out;

layout(std430, binding = 1) buffer boatPositionsSSBO{
    vec4 boatPositions[];
};

layout(std430, binding = 6) buffer boatAngularPositionsSSBO {
    vec4 boatAngularPositions[];
};

uniform int boatIndex;
uniform mat4 m_mvp;
uniform float waveHeight;
uniform mat4 m_vp;
uniform mat4 m_model;
uniform sampler2D oceanHeightmap;

out vec4 colorIn;

// Returns a texture coordinate by normalizing a 2D vector
// within x and z [-32, 32]
vec2 heightmapCoordinate(float x, float z) {
    return vec2((x + 32)/64, (z + 32)/64);
}

float waveHeightAtPoint(float x, float z) {
    return texture(oceanHeightmap, heightmapCoordinate(x, z)).r * waveHeight;
}

// Returns a vector that is normal to 3 points
//vec3 normal(vec3 A, vec3 B, vec3 C) {
//    vec3 U = B - A;
//    vec3 V = C - A;
//    vec3 N = cross(U, V);
//
//    if (N != vec3(0.0)) {
//        N = normalize(N);
//    }
//    return N;
//}

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
// !!!Really, I have no idea why this original commented method is half as complex as it is here.
// Newer version seems to work just fine and is infinitely simpler.!!!
//float surfaceDistance(vec3 vertex) {
//    vec4 pl = oceanPlane(vertex.x, vertex.z);
//
//    float oceanAltitude = -((vertex.x * pl.x) + (vertex.z * pl.z) + pl.w) / pl.y;
//    return vertex.y - oceanAltitude;
//}

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

// Apply boat position and rotation transforms to a vector
vec3 applyBoatTransforms(vec3 vertex) {
    return rotate(vertex, boatAngularPositions[boatIndex].xyz) + boatPositions[boatIndex].xyz;
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

void main() {
    vec3 A = (m_model * gl_in[0].gl_Position).xyz;
    vec3 B = (m_model * gl_in[1].gl_Position).xyz;
    vec3 C = (m_model * gl_in[2].gl_Position).xyz;
    A = applyBoatTransforms(A);
    B = applyBoatTransforms(B);
    C = applyBoatTransforms(C);

    float ASurfaceDist = surfaceDistance(A);
    float BSurfaceDist = surfaceDistance(B);
    float CSurfaceDist = surfaceDistance(C);

    bool AIsSubmerged = ASurfaceDist < 0;
    bool BIsSubmerged = BSurfaceDist < 0;
    bool CIsSubmerged = CSurfaceDist < 0;


    int totalSubmergedVertices = (AIsSubmerged ? 1 : 0) + (BIsSubmerged ? 1 : 0) + (CIsSubmerged ? 1 : 0);
    int workingTriangles = 1;
    vec3[9] vertices;

    vec3[3] baseVerticesAux = {A, B, C};
    float[3] baseDistsAux = {ASurfaceDist, BSurfaceDist, CSurfaceDist};

    switch (totalSubmergedVertices) {
        case 0: {
            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;
        }
            break;
        case 1: {
            workingTriangles = 3;

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
            workingTriangles = 3;
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
            vertices[0] = A;
            vertices[1] = B;
            vertices[2] = C;
        }
            break;
        default: break;
    }
 
    transformAndEmitVertices(vertices, workingTriangles);
    //vertices[0] = A;
    //vertices[1] = B;
    //vertices[2] = C;
    //transformAndEmitVertices(vertices, 1);
}