#version 450

layout (std430, binding = 2) buffer forcesSSBO {
    vec4 forces[];
};

out vec4 color;

void main() {
    //color = colorIn;
    //color = vec4(forces[0].xyz, 1.0);
    color = vec4(1.0);
}