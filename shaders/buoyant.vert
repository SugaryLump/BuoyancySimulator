#version 450

in vec3 vertex;

uniform mat4 m_mvp;

void main() {
    gl_Position = vec4(vertex, 1.0);
}