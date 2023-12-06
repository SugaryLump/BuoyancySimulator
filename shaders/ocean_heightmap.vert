#version 450
#define PI 3.1415926538

in vec3 vertex;

out vec3 frag_vertex;

void main() {
    frag_vertex = vertex;
    gl_Position = vec4(vertex, 1.0);
}