#include <GL/glew.h>
#include <GL/freeglut.h>

#include "main.hpp"
#include "camera.hpp"
#include "models.hpp"
#include "shader.hpp"

#include <iostream>
#include <vec4.hpp>
#include <memory>
#include <tiny_obj_loader.h>
#include <mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

using namespace std;
using namespace glm;

#define WAVE_HEIGHT 0.3f;

shared_ptr<Camera> camera;
shared_ptr<Shader> buoyantShader;

auto buoyantModels = vector<shared_ptr<Model>>();

GLuint oceanHeightmapFramebuffer;
GLuint oceanHeightmapTexture;
shared_ptr<Shader> oceanHeightmapShader;
shared_ptr<Model> oceanGrid;

shared_ptr<Shader> oceanShader;
shared_ptr<Model> oceanPlane;

GLuint verticesSSBO, positionsSSBO, angularPositionsSSBO;

GLuint elapsedTime = 0;
GLuint frameFrequency = 0;

void updateEllapsedTime() {
    GLuint newTime = glutGet(GLUT_ELAPSED_TIME);
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

    // Enable wireframe mode (debugging purposes, remove later)
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // Draw all buoyant models with the buoyant shader 
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

    buoyantShader->BindShader();

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, oceanHeightmapTexture);
    {
        mat4 m_vp = camera->GetViewProjectionMatrix();
        GLuint m_vpLocation = glGetUniformLocation(buoyantShader->shaderProgram, "m_vp");
        glUniformMatrix4fv(m_vpLocation, 1, GL_FALSE, value_ptr(m_vp));

        GLfloat waveHeight = WAVE_HEIGHT;
        GLuint waveHeightLocation = glGetUniformLocation(buoyantShader->shaderProgram, "waveHeight");
        glUniform1f(waveHeightLocation, waveHeight);
        for (int i = 0; i < buoyantModels.size(); i++) {
            mat4 m_model = buoyantModels[i]->GetModelMatrix();
            GLuint m_modelLocation = glGetUniformLocation(buoyantShader->shaderProgram, "m_model");
            glUniformMatrix4fv(m_modelLocation, 1, GL_FALSE, value_ptr(m_model));

            mat4 m_mvp = m_vp * m_model;
            GLuint m_mvpLocation = glGetUniformLocation(buoyantShader->shaderProgram, "m_mvp");
            glUniformMatrix4fv(m_mvpLocation, 1, GL_FALSE, value_ptr(m_mvp));

            GLuint boatIndex = i;
            GLuint boatIndexLocation = glGetUniformLocation(buoyantShader->shaderProgram, "boatIndex");
            glUniform1ui(boatIndexLocation, (boatIndex));

            GLuint texLocation = glGetUniformLocation(buoyantShader->shaderProgram, "oceanHeightmap");
            glUniform1i(texLocation, GL_TEXTURE0);

            buoyantModels[i]->draw();
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
    float positions_aux[buoyantModels.size() * 4];
    float angularPositions_aux[buoyantModels.size() * 4];

    for (int i = 0; i < buoyantModels.size(); i++) {
        positions_aux[i * 4] = 0.0f;
        positions_aux[i * 4 + 1] = 0.0f;
        positions_aux[i * 4 + 2] = 0.0f;
        positions_aux[i * 4 + 3] = 0.0f;

        angularPositions_aux[i * 4] = 0.0f;
        angularPositions_aux[i * 4 + 1] = 0.0f;
        angularPositions_aux[i * 4 + 2] = 0.0f;
        angularPositions_aux[i * 4 + 3] = 0.0f;
    }

    //glGenBuffers(1, &verticesSSBO);
    //glBindBuffer(GL_SHADER_STORAGE_BUFFER, verticesSSBO);
    //glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), NULL, GL_DYNAMIC_DRAW);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, verticesSSBO);
    
    glGenBuffers(1, &positionsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, positionsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), positions_aux, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, positionsSSBO);

    glGenBuffers(1, &angularPositionsSSBO);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, positionsSSBO);
    glBufferData(GL_SHADER_STORAGE_BUFFER, buoyantModels.size() * sizeof(vec4), angularPositions_aux, GL_DYNAMIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 6, angularPositionsSSBO);

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
    buoyantModels.push_back(make_shared<Model>(argv[1]));

    // Load shaders
    oceanShader = make_shared<Shader>("shaders/ocean");
    oceanHeightmapShader = make_shared<Shader>("shaders/ocean_heightmap");
    buoyantShader = make_shared<Shader>("shaders/buoyant");

    // Initialize SSBOs
    initSSBOs();

    // Begin main loop
    cout << "Entering main render loop\n" << flush;
    glutMainLoop();

    return 0;
}