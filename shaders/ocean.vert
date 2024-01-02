#version 450

in vec3 vertex;
in vec2 tex_coords;


uniform mat4 m_model;
uniform mat4 m_vp;
uniform sampler2D height_map;
uniform float waveHeight;

void main() {
    vec4 translated = m_model * vec4(vertex, 1.0);
    translated.y = texture(height_map, tex_coords).r * waveHeight;
    gl_Position = m_vp * translated;
}