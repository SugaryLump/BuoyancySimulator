#version 450
#define PI 3.1415926538

in vec3 vertex;
out float height;

uniform int ellapsedTime;

void main() {
    float x = (vertex.x + 1) / 2 * PI;
    height = sin(vertex.x + ellapsedTime * (PI / 2000));
    gl_Position = vec4(vertex, 1.0);
}