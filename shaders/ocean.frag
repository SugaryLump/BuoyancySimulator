#version 450

in vec2 tex_coords_frag;

out vec4 color;

uniform sampler2D height_map;


void main() {
    vec4 text = texture(height_map, vec2(0.5, 0.5));
    color = vec4(text.rgb, 1.0);
}