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

shared_ptr<Camera> camera;
shared_ptr<Shader> cameraShader;

auto buoyantModels = vector<shared_ptr<Model>>();

GLuint oceanHeightmapFramebuffer;
shared_ptr<Shader> oceanHeightmapShader;
shared_ptr<Shader> oceanShader;
shared_ptr<Model> oceanPlane;

void render() {
    // Update camera
    camera->UpdatePosition();

    // Set background colour
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    // Draw ocean height map
    oceanHeightmapShader->BindShader();
    glBindFramebuffer(GL_FRAMEBUFFER, oceanHeightmapFramebuffer);
    glViewport(0, 0, FRAMEBUFFER_SIZE, FRAMEBUFFER_SIZE);
    oceanPlane->draw();
    

    // Draw all models with the camera shader
    glViewport(0, 0, camera->windowWidth, camera->windowHeight);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);    
    cameraShader->BindShader();
    for (auto model : buoyantModels) {
        mat4 mvp = camera->GetViewProjectionMatrix() * model->GetModelMatrix();
        GLuint matLocation = glGetUniformLocation(cameraShader->shaderProgram, "m_mvp");
        glUniformMatrix4fv(matLocation, 1, GL_FALSE, value_ptr(mvp));
        model->draw();
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

GLuint initTextureRenderTarget() {
    GLuint frameBuffer = 0;
    glGenFramebuffers(1, &frameBuffer);
    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, FRAMEBUFFER_SIZE, FRAMEBUFFER_SIZE, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture, 0);

    GLenum DrawBuffers[1] = {GL_COLOR_ATTACHMENT0};
    glDrawBuffers(1, DrawBuffers);

    if(glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        cout << "Error initializing frame buffer." << endl << flush;
        exit(-1);
    }
    
    return frameBuffer;
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

    // Face culling
    glEnable(GL_CULL_FACE);

    // Initialize oceanHeightmap render target
    oceanHeightmapFramebuffer = initTextureRenderTarget();

    // Load models
    oceanPlane = make_shared<Model>("models/quad.obj");
    buoyantModels.push_back(make_shared<Model>(argv[1]));

    // Load shaders
    oceanHeightmapShader = make_shared<Shader>("shaders/ocean_heightmap");
    cameraShader = make_shared<Shader>("shaders/camera");

    // Begin main loop
    cout << "Entering main render loop\n" << flush;
    glutMainLoop();

    return 0;
}