#version 450

in vec4 colorIn;
out vec4 color;

void main() {
    //color = vec4(1.0, 0.0, 0.0, 1.0);
    color = colorIn;
}