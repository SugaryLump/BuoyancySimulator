#version 450

in float height;

layout(location = 0) out vec4 color;

void main() {
    color = vec4(height, 1.0, 1.0, 1.0);
}