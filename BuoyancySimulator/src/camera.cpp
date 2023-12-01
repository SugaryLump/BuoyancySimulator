#include "camera.hpp"
#include <GL/freeglut.h>
#include <vec4.hpp>

using namespace glm;

Camera::Camera() {
  this->near = 0.1;
  this->far = 100;
  this->fov = 90;

  this->position = vec4(0, 0, 0, 0);
  this->lookAt = vec4(0, 1, 0, 0);
  this->up = vec4(0, 0, 1, 0);
}

void Camera::UpdateMatrices(float ratio) {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //gluPerspective(this->fov, 1, this->near, this->far);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  //gluLookAt(this->position.x, this->position.y, this->position.z,
  //          this->lookAt.x, this->lookAt.y, this->lookAt.z,
  //          this->up.x, this->up.y, this->up.z);
}