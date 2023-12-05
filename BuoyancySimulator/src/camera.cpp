#include <GL/glew.h>
#include <GL/freeglut.h>

#include "camera.hpp"

#include <vec4.hpp>
#include <mat4x4.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/transform.hpp>
#include <math.h>
#include <cstring>
#include <string>
#include <iostream>

using namespace glm;
using namespace std;

Camera::Camera() {
  this->aspectRatio = 1;
  this->windowWidth = 500;
  this->windowHeight = this->windowWidth * this->aspectRatio;

  this->near = 0.001;
  this->far = 1000;
  this->fov = M_PI_4;
  this->position = {0, 0, 0};
  this->pitch = 0;
  this->yaw = 0;
  memset(this->keyboardState, 0, sizeof(bool) * 256);

  this->UpdateVectorsAndMatrices();
}

void Camera::Update(GLuint frameFrequency) {
  this->UpdatePosition(frameFrequency);
  this->UpdateVectorsAndMatrices();
}

void Camera::UpdateVectorsAndMatrices() {
  this->lookVector = {
    sin(this->yaw) * cos(this->pitch),
    sin(this->pitch),
    -cos(this->yaw) * cos(this->pitch)
  };
  this->rightVector = {
    sin(this->yaw + M_PI_2) * cos(this->pitch),
    0,
    -cos(this->yaw + M_PI_2) * cos(this->pitch)
  };
  this->upVector = {
    sin(this->yaw) * cos(this->pitch + M_PI_2),
    sin(this->pitch + M_PI_2),
    -cos(this->yaw) * cos(this->pitch + M_PI_2)
  };

  this->viewMatrix = lookAt(
    this->position,
    this->lookVector + this->position,
    this->upVector
  );

  this->projectionMatrix = perspective(
    this->fov,
    this->aspectRatio,
    this->near,
    this->far
  );
  
  this->viewProjectionMatrix = this->projectionMatrix * this->viewMatrix;
}

void Camera::PressKey(unsigned char key) {
  this->keyboardState[key] = true;
}

void Camera::LiftKey(unsigned char key) {
  this->keyboardState[key] = false;
}

void Camera::UpdatePosition(GLuint frameFrequency) {
  vec3 translation = {0, 0, 0};

  if (this->keyboardState['w'] || this->keyboardState['W']) {
    translation += this->lookVector;
  }
  if (this->keyboardState['a'] || this->keyboardState['A']) {
    translation -= this->rightVector;
  }
  if (this->keyboardState['s'] || this->keyboardState['S']) {
    translation -= this->lookVector;
  }
  if (this->keyboardState['d'] || this->keyboardState['D']) {
    translation += this->rightVector;
  }
  if (this->keyboardState[' ']) {
    translation += vec3{0, 1, 0};
  }
  if (this->keyboardState['v'] || this->keyboardState['V']) {
    translation -= vec3{0, 1, 0};
  }

  if (translation != vec3{0, 0, 0}) {
    this->position += normalize(translation) * MOVEMENT_SPEED * (frameFrequency/1000.0f);
    this->UpdateVectorsAndMatrices();
  }
}

// x and y are offset from center
void Camera::ProcessMouseMotion(float x, float y) {
  this->yaw += x * CAMERA_SENSITIVITY;
  this->pitch -= y * CAMERA_SENSITIVITY;

  if (this->yaw < 0) {
    this->yaw += 2 * M_PI;
  }
  else if (this->yaw > 2 * M_PI) {
    this->yaw -= 2 * M_PI;
  }

  if (this->pitch > M_PI_2 - 0.01) {
    this->pitch = M_PI_2 - 0.01;
  }
  else if (this->pitch < -M_PI_2 + 0.01) {
    this->pitch = -M_PI_2 + 0.01;
  }
}

mat4 Camera::GetViewMatrix() {
  return this->viewMatrix;
}

mat4 Camera::GetProjectionMatrix() {
  return this->projectionMatrix;
}

mat4 Camera::GetViewProjectionMatrix() {
  return this->viewProjectionMatrix;
}