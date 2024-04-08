#version 450

layout (std430, binding = 2) buffer forcesSSBO {
    vec4 forces[];
};

out vec4 color;

void main() {
    //color = colorIn;
    //color = vec4(forces[0].xyz, 1.0);
    //color = vec4(1.0, forces[0].y, forces[0].y, 1.0);
    color = vec4(vec3(1.0), 1);
}