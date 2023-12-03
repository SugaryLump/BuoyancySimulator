#version 450

in vec3 vertex;

uniform mat4 m_mvp;

void main() {
    gl_Position = m_mvp * vec4(vertex, 1);
}