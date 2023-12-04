#version 450

in vec3 vertex;
out vec3 world_pos;
void main() {
    world_pos = vertex;
    gl_Position = vec4(vertex, 1.0);
}