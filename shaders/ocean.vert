#version 450

in vec3 vertex;
in vec2 tex_coords;

out vec2 tex_coords_frag;
uniform mat4 m_mvp;
uniform sampler2D height_map;

void main() {
    vec4 translated = vec4(vertex.x, texture(height_map, tex_coords).r, vertex.z, 1.0);
    tex_coords_frag = tex_coords;
    gl_Position = m_mvp * vec4(vertex, 1.0);
    //gl_Position = m_mvp * translated;
}