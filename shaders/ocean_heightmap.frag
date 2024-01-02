#version 450
#define PI 3.1415926538

in vec3 frag_vertex;

uniform uint elapsedTime;

layout(location = 0) out vec4 color;

void main() {
    float x = (frag_vertex.x + 1) * 4 * PI;
    float height = (sin(x + elapsedTime * (PI / 3000)) + 1) / 2;
    color = vec4(height, 0, 0, 1.0);
}