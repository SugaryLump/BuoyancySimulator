#version 450

in vec3 world_pos;

layout(location = 0) out vec4 color;

void main() {
    color = vec4(world_pos.x, world_pos.y, world_pos.z, 1.0);
}