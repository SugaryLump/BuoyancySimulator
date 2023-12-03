#include <GL/glew.h>
#include <GL/freeglut.h>

#include "main.hpp"
#include "camera.hpp"
#include "models.hpp"

#include <iostream>
#include <vec4.hpp>
#include <memory>
#include <tiny_obj_loader.h>

using namespace std;
using namespace glm;

float ratio = 1;
auto cam = make_shared<Camera>();
auto models = vector<shared_ptr<Model>>();

void testRender() {
    /*
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-0.5f, -0.5f, 1.0f);
    glVertex3f(0.5f, -0.5f, 1.0f);
    glVertex3f(0.5f, 0.5f, 1.0f);
    glVertex3f(-0.5f, 0.5f, 1.0f);
    */
    //glBindVertexArray(*vao);
    //glDrawRangeElements(GL_TRIANGLES, 0, 3, 6, GL_UNSIGNED_SHORT, NULL);
}

void render() {
    cam->UpdateMatrices(ratio);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    for (auto model : models) {
        model->draw();
    }

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

    // Init glew
    glewInit();

    // Load model
    models.push_back(make_shared<Model>(argv[1]));

    // Begin main loop
    cout << "Entering main render loop\n" << flush;
    glutMainLoop();

    return 0;
}