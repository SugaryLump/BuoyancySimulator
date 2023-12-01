#include "main.hpp"
#include "camera.hpp"

#include <iostream>
#include <GL/freeglut.h>
#include <vec4.hpp>
#include <memory>

using namespace std;
using namespace glm;

float ratio = 1;
shared_ptr<Camera> cam = make_shared<Camera>();

void testRender() {
    glBegin(GL_QUADS);
    glColor3f(1.0f, 0.0f, 1.0f);
    glVertex3f(-0.5f, -0.5f, 1.0f);
    glVertex3f(0.5f, -0.5f, 1.0f);
    glVertex3f(0.5f, 0.5f, 1.0f);
    glVertex3f(-0.5f, 0.5f, 1.0f);
    glEnd();
}

void render() {
    cam->UpdateMatrices(ratio);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    testRender();

    glutSwapBuffers();
}

void resizeWindow(int x, int y) {
    if (y == 0) {
        y = 1;
    }

    glViewport(0, 0, x, y);

    ratio = (float)x / y;
}

int main(int argc, char* argv[]) {
    // GLUT initialization and window settings
    cout << "Initializing GLUT\n" << flush;
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA | GLUT_MULTISAMPLE);

    glutInitWindowPosition(100, 100);
    glutInitWindowSize(512, 512);
    glutCreateWindow("Buoyancy Simulator");

    // Callbacks
    glutDisplayFunc(render);
    glutReshapeFunc(resizeWindow);

    cout << "Entering main render loop\n" << flush;
    glutMainLoop();
    return 0;
}