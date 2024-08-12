#include <GL/glew.h>
#include <GL/freeglut.h>

#include <iostream>
#include <vec4.hpp>
#include <memory>
#include <tiny_obj_loader.h>
#include <mat3x3.hpp>
#include <mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#define TINYOBJLOADER_IMPLEMENTATION
#define GLM_ENABLE_EXPERIMENTAL

#include "main.hpp"
#include "camera.hpp"
#include "models.hpp"
#include "shader.hpp"
#include "buoyant.hpp"

using namespace std;
using namespace glm;

#define WAVE_HEIGHT 0.3f;

shared_ptr<Camera> camera;

GLuint oceanHeightmapFramebuffer;
GLuint oceanHeightmapTexture;
shared_ptr<Shader> oceanHeightmapShader;
shared_ptr<Model> oceanGrid;

shared_ptr<Shader> oceanShader;
shared_ptr<Model> oceanPlane;

GLuint boatPositionsSSBO, forcesSSBO, boatVelocitiesSSBO,
       boatAngularVelocitiesSSBO, boatAngularPositionsSSBO, torquesSSBO,
       boatOldTriangleVelocitiesSSBO, boatOldTriangleAngularVelocitiesSSBO, boatOldSubmergedAreasSSBO;
shared_ptr<Shader> buoyantCalcsShader;
auto buoyantModels = vector<shared_ptr<Buoyant>>();
#define REDUCE_SHADER_GROUP_SIZE 16
size_t maxTriangleCount;
shared_ptr<Shader> forcesReductionShader;
shared_ptr<Shader> torquesReductionShader;
shared_ptr<Shader> buoyantApplicationShader;
shared_ptr<Shader> forceVisualizationShader;
shared_ptr<Shader> cameraShader;

GLuint elapsedTime = 0;
GLuint frameFrequency = 0;
GLuint oldFrameFrequency = 0;

void updateEllapsedTime() {
    GLuint newTime = glutGet(GLUT_ELAPSED_TIME);
    oldFrameFrequency = frameFrequency;
    frameFrequency = newTime - elapsedTime;
    elapsedTime = newTime;
}

void render() {
    // Get frame frequency
    updateEllapsedTime();

    // Update camera
    camera->Update(frameFrequency);

    // Set background colour and clear frame
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Enable fill mode
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // Calculate ocean height map
    glDisable(GL_DEPTH_TEST);

    oceanHeightmapShader->BindShader();

    glBindFramebuffer(GL_FRAMEBUFFER, oceanHeightmapFramebuffer);
    glViewport(0, 0, FRAMEBUFFER_SIZE, FRAMEBUFFER_SIZE);

    {
        GLuint elapsedTimeLocation = glGetUniformLocation(oceanHeightmapShader->shaderProgram, "elapsedTime");
        glUniform1ui(elapsedTimeLocation, elapsedTime);

        oceanGrid->draw();
    }

    // Set framebuffer back to default for drawing on scren
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, camera->windowWidth, camera->windowHeight);

    // Enable wireframe mode (debugging purposes; remove later)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Run the buoyant shaders on buoyant models 
    {
        // Get shader uniform locations and calculate needed data
        mat4 m_vp = camera->GetViewProjectionMatrix();
        GLuint m_vpLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "m_vp");
        GLuint m_scaleRotationLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "m_scale_rotation");
        GLuint m_translationLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "m_translation");
        GLfloat waveHeight = WAVE_HEIGHT;
        GLuint waveHeightLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "waveHeight");
        GLuint boatMassLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "boatMass");
        GLuint boatTotalAreaLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "boatTotalArea");
        GLuint maxTrianglesLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "maxTriangles");
        GLuint boatIndexLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "boatIndex");
        GLuint boatLengthLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "boatLength");
        GLuint boatCenterOfMassCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "boatCenterOfMass");
        GLuint texLocationCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "oceanHeightmap");
        GLuint oldDeltaTimeCalcs = glGetUniformLocation(buoyantCalcsShader->shaderProgram, "oldDeltaTime");


        GLuint m_scaleRotationLocationVis = glGetUniformLocation(forceVisualizationShader->shaderProgram, "m_scale_rotation");
        GLuint m_translationLocationVis = glGetUniformLocation(forceVisualizationShader->shaderProgram, "m_translation");
        GLuint m_vpLocationVis = glGetUniformLocation(forceVisualizationShader->shaderProgram, "m_vp");
        GLuint boatIndexLocationVis = glGetUniformLocation(forceVisualizationShader->shaderProgram, "boatIndex");


        GLuint totalForcesLocationReducs = glGetUniformLocation(forcesReductionShader->shaderProgram, "totalForces");
        GLuint totalTorquesLocationReducs = glGetUniformLocation(torquesReductionShader->shaderProgram, "totalTorques");


        GLuint boatIndexLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "boatIndex");
        GLuint boatMassLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "boatMass");
        GLuint boatInertiaModifierLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "boatInertiaModifier");
        GLuint deltaTimeLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "deltaTime");
        GLuint m_vpLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "m_vp");
        GLuint m_scaleRotationLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "m_scale_rotation");
        GLuint m_translationLocation = glGetUniformLocation(buoyantApplicationShader->shaderProgram, "m_translation");
        for (int i = 0; i < buoyantModels.size(); i++) {
            // Run hydrodynamic/hydrostatic force calculations
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);

            buoyantCalcsShader->BindShader();

            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, oceanHeightmapTexture);

            glUniformMatrix4fv(m_vpLocationCalcs, 1, GL_FALSE, value_ptr(m_vp));
            
            glUniform1f(waveHeightLocationCalcs, waveHeight);

            mat4 m_scale_rotation = buoyantModels[i]->GetModel().GetScaleRotationMatrix();
            glUniformMatrix4fv(m_scaleRotationLocationCalcs, 1, GL_FALSE, value_ptr(m_scale_rotation));

            mat4 m_translation = buoyantModels[i]->GetModel().GetTranslationMatrix();
            glUniformMatrix4fv(m_translationLocationCalcs, 1, GL_FALSE, value_ptr(m_translation));

            GLfloat boatMass = buoyantModels[i]->GetMass();
            glUniform1f(boatMassLocationCalcs, boatMass);

            GLfloat totalArea = buoyantModels[i]->GetModel().GetTotalArea();
            glUniform1f(boatTotalAreaLocationCalcs, totalArea);

            GLfloat boatLength = buoyantModels[i]->GetModel().GetLength();
            glUniform1f(boatLengthLocationCalcs, boatLength);

            vec3 boatCenterOfMass = buoyantModels[i]->GetCenterOfMass();
            glUniform3f(boatCenterOfMassCalcs, boatCenterOfMass.x, boatCenterOfMass.y, boatCenterOfMass.z);

            glUniform1ui(oldDeltaTimeCalcs, oldFrameFrequency);

            GLuint maxTriangles = maxTriangleCount;
            glUniform1ui(maxTrianglesLocationCalcs, maxTriangles);

            GLuint boatIndex = i;
            glUniform1ui(boatIndexLocationCalcs, (boatIndex));

            glUniform1i(texLocationCalcs, GL_TEXTURE0);

            buoyantModels[i]->GetModel().draw();

            // Run voxel debug visualization
            {
                cameraShader->BindShader();
                
                mat4 m_mvpdebug = camera->GetViewProjectionMatrix();
                GLuint m_mvpdebuglocation = glGetUniformLocation(cameraShader->shaderProgram, "m_mvp");
                glUniformMatrix4fv(m_mvpdebuglocation, 1, GL_FALSE, value_ptr(m_mvpdebug));
                
                buoyantModels[i]->GetModel().drawVoxelsDebug();
            }

            // Run force visualization shader
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_CULL_FACE);

            forceVisualizationShader->BindShader();

            glUniformMatrix4fv(m_vpLocationVis, 1, GL_FALSE, value_ptr(m_vp));
            glUniformMatrix4fv(m_scaleRotationLocationVis, 1, GL_FALSE, value_ptr(m_scale_rotation));
            glUniformMatrix4fv(m_translationLocationVis, 1, GL_FALSE, value_ptr(m_translation));

            glUniform1ui(boatIndexLocationVis, boatIndex);

            buoyantModels[i]->GetModel().draw();

            // Run force and torque SSBO sum reduction
            // The "q = (dividend + divisor - 1) / dividor" style divisions
            // here are just a hacky way of getting the ceiling of fractional
            // when we're only working with integers
            for (GLuint totalElements = buoyantModels[i]->GetModel().GetTriangleCount() * 3; totalElements > 1;
                 totalElements = (totalElements + REDUCE_SHADER_GROUP_SIZE * 2 - 1) / (REDUCE_SHADER_GROUP_SIZE * 2)) {
                GLuint totalWorkGroups = (totalElements + REDUCE_SHADER_GROUP_SIZE * 2 - 1) / (REDUCE_SHADER_GROUP_SIZE * 2);
                
                forcesReductionShader->BindShader();
                glUniform1ui(totalForcesLocationReducs, totalElements);
                glDispatchCompute(totalWorkGroups, 1, 1);
                
                torquesReductionShader->BindShader();
                
                glUniform1ui(totalTorquesLocationReducs, totalElements);
                glDispatchCompute(totalWorkGroups, 1, 1);
                
                glMemoryBarrier(GL_ALL_BARRIER_BITS);
            }

            // Render the models with the force application shader
            glEnable(GL_DEPTH_TEST);
            glEnable(GL_CULL_FACE);

            buoyantApplicationShader->BindShader();


            glUniform1i(boatIndexLocation, i);
            glUniform1f(boatMassLocation, buoyantModels[i]->GetMass());
            glUniform1f(boatInertiaModifierLocation, buoyantModels[i]->GetInertiaModifier());
            glUniform1ui(deltaTimeLocation, frameFrequency);
            glUniformMatrix4fv(m_vpLocation, 1, GL_FALSE, value_ptr(m_vp));
            glUniformMatrix4fv(m_scaleRotationLocation, 1, GL_FALSE, value_ptr(m_scale_rotation));
            glUniformMatrix4fv(m_translationLocation, 1, GL_FALSE, value_ptr(m_translation));
                
            buoyantModels[i]->GetModel().draw();
            glMemoryBarrier(GL_ALL_BARRIER_BITS);
        }
    }
    

    // Draw ocean
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);

    oceanShader->BindShader();
    
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, oceanHeightmapTexture);
    {
        mat4 m_vp = camera->GetViewProjectionMatrix();
        GLuint m_vpLocation = glGetUniformLocation(oceanShader->shaderProgram, "m_vp");
        glUniformMatrix4fv(m_vpLocation, 1, GL_FALSE, value_ptr(m_vp));

        mat4 m_model = oceanPlane->GetModelMatrix();
        GLuint m_modelLocation = glGetUniformLocation(oceanShader->shaderProgram, "m_model");
        glUniformMatrix4fv(m_modelLocation, 1, GL_FALSE, value_ptr(m_model));

        GLuint texLocation = glGetUniformLocation(oceanShader->shaderProgram, "height_map");
        glUniform1i(texLocation, GL_TEXTURE0);

        GLfloat waveHeight = WAVE_HEIGHT;
        GLuint waveHeightLocation = glGetUniformLocation(oceanShader->shaderProgram, "waveHeight");
        glUniform1f(waveHeightLocation, waveHeight);

        oceanPlane->draw();
    }

    // Force redraw
    glutPostRedisplay();
    glutSwapBuffers();
}

void resizeWindow(int x, int y) {
    if (y == 0) {
        y = 1;
    }

    camera->windowWidth = x;
    camera->windowHeight = y;
    camera->aspectRatio = (float)x/y;
}

void pressKey(unsigned char key, int x, int y) {
    camera->PressKey(key);
    glutWarpPointer(camera->windowWidth/2, camera->windowHeight/2);
}

void liftKey(unsigned char key, int x, int y) {
    camera->LiftKey(key);
}

void moveMouse(int x, int y) {
    x = x - camera->windowWidth/2;
    y = y - camera->windowHeight/2;
    if (x != 0 || y != 0) {
        camera->ProcessMouseMotion(x, y);
        glutSetCursor(GLUT_CURSOR_NONE);
        glutWarpPointer(camera->windowWidth/2, camera->windowHeight/2);
    }
}

void initOceanHeightmapRenderTarget() {
    oceanHeightmapFramebuffer = 0;
    glGenFramebuffers(1, &oceanHeightmapFramebuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, oceanHeightmapFramebuffer);

    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &oceanHeightmapTexture);
    glBindTexture(GL_TEXTURE_2D, oceanHeightmapTexture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, FRAMEBUFFER_SIZE, FRAMEBUFFER_SIZE, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, oceanHeightmapTexture, 0);

    GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, DrawBuffers);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        cout << "Error initializing frame buffer." << endl << flush;
        exit(-1);
    }
}

void initSSBOs() {
    maxTriangleCount = buoyantModels[0]->GetModel().GetTriangleCount();
    for (int i = 1; i < buoyantModels.size(); i++) {
        size_t triangleCount = buoyantModels[i]->GetModel().GetTriangleCount();
        if (triangleCount > maxTriangleCount) {
            maxTriangleCount = triangleCount;
        }
    }
    float nBoatsZeroVec4s[buoyantModels.size() * 4];
    memset(nBoatsZeroVec4s, 0, sizeof(nBoatsZeroVec4s));
    float maxSplitTrianglesZeroVec4s[maxTriangleCount * 3 * 4];
    memset(maxSplitTrianglesZeroVec4s, 0, sizeof(maxSplitTrianglesZeroVec4s));
    float maxTrianglesZeroVec4s[maxTriangleCount * 4];
    memset(maxTrianglesZeroVec4s, 0, sizeof(maxTrianglesZeroVec4s));
    float maxTrianglesZeroFloats[maxTriangleCount];
    memset(maxTrianglesZeroFloats, 0, sizeof(maxTrianglesZeroFloats));
    float maxTrianglesZeroFloatsPerBoat[buoyantModels.size() * maxTriangleCount];
    memset(maxTrianglesZeroFloatsPerBoat, 0, sizeof(maxTrianglesZeroFloatsPerBoat));
    float maxTrianglesZeroVec4sPerBoat[buoyantModels.size() * maxTriangleCount * 4];
    memset(maxTrianglesZeroVec4sPerBoat, 0, sizeof(maxTrianglesZeroVec4sPerBoat));

    // GLuint boatPositionsSSBO, forcesSSBO, boatVelocitiesSSBO,
    //   boatAngularPositionsSSBO, boatAngularPositionsSSBO, torquesSSBO;

    glGenBuffers(1, &boatPositionsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatPositionsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), nBoatsZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, boatPositionsSSBO);

    glGenBuffers(1, &forcesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, forcesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, maxTriangleCount * sizeof(vec4) * 3, maxSplitTrianglesZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, forcesSSBO);

    glGenBuffers(1, &boatVelocitiesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatVelocitiesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), nBoatsZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, boatVelocitiesSSBO);

    glGenBuffers(1, &boatAngularVelocitiesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatAngularVelocitiesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), nBoatsZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, boatAngularVelocitiesSSBO);

    glGenBuffers(1, &boatOldTriangleVelocitiesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatOldTriangleVelocitiesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * maxTriangleCount * sizeof(vec4), nBoatsZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 9, boatOldTriangleVelocitiesSSBO);

    glGenBuffers(1, &boatOldSubmergedAreasSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatOldSubmergedAreasSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * maxTriangleCount * sizeof(float), maxTrianglesZeroFloats, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 10, boatOldSubmergedAreasSSBO);

    glGenBuffers(1, &boatAngularPositionsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, boatAngularPositionsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), nBoatsZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, boatAngularPositionsSSBO);

    glGenBuffers(1, &torquesSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, torquesSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, maxTriangleCount * sizeof(vec4) * 3, maxSplitTrianglesZeroVec4s, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 8, torquesSSBO);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
}

int main(int argc, char* argv[]) {
    // GLUT initialization and context creation
    cout << "Initializing GLUT\n" << flush;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE);
    glutInitContextVersion(4,6);
    glutInitContextProfile(GLUT_CORE_PROFILE);
    glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);

    // Create camera and window
    camera = make_shared<Camera>();
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(camera->windowWidth, camera->windowHeight);
    glutCreateWindow("Buoyancy Simulator");

    // Callbacks
    glutDisplayFunc(render);
    glutReshapeFunc(resizeWindow);
    glutKeyboardFunc(pressKey);
    glutKeyboardUpFunc(liftKey);
    glutPassiveMotionFunc(moveMouse);

    // Init glew
    glewInit();

    // Alpha blending
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );

    // Initialize oceanHeightmap render target
    initOceanHeightmapRenderTarget();

    // Load models
    oceanGrid = make_shared<Model>("models/ocean_grid.obj");
    oceanPlane = make_shared<Model>("models/ocean_plane.obj");
    oceanPlane->SetScale(vec3(64, 0, 64));
    if (argc > 1) {
        for (int i = 1; i > argc; i++) {
            buoyantModels.push_back(make_shared<Buoyant>(argv[i], true));
        }
    }
    else {
        buoyantModels.push_back(make_shared<Buoyant>());
    }

    // Load shaders
    oceanShader = make_shared<Shader>("shaders/ocean");
    oceanHeightmapShader = make_shared<Shader>("shaders/ocean_heightmap");
    buoyantCalcsShader = make_shared<Shader>("shaders/buoyant_calcs");
    forcesReductionShader = make_shared<Shader>("shaders/forces_reduction");
    torquesReductionShader = make_shared<Shader>("shaders/torques_reduction");
    buoyantApplicationShader = make_shared<Shader>("shaders/buoyant_application");
    forceVisualizationShader = make_shared<Shader>("shaders/force_visualization");
    cameraShader = make_shared<Shader>("shaders/camera");

    // Initialize SSBOs
    initSSBOs();

    // Begin main loop
    cout << "Entering main render loop\n" << flush;
    glutMainLoop();

    return 0;
}