#include <GL/glew.h>
#include <GL/freeglut.h>

#include "shader.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>
#include <vector>

using namespace std;

vector<string> shaderTypes = {
    "vert",
    "geom",
    "frag",
    "comp"
};
vector<int> glShaderTypes = {
    GL_VERTEX_SHADER,
    GL_GEOMETRY_SHADER,
    GL_FRAGMENT_SHADER,
    GL_COMPUTE_SHADER,
};

string readFile(string fileName)
{
    ifstream file = ifstream(fileName);
    if (!file.is_open())
    {
        return "";
    }
    stringstream strStream;
    strStream << file.rdbuf();
    file.close();

    return strStream.str();
}

Shader::Shader(string shaderName)
{
    this->shaderProgram = 0;
    vector<GLuint> compiledShaders;

    for (int i = 0; i < shaderTypes.size(); i++)
    {
        string shaderType = shaderTypes[i];
        string shaderString = readFile(shaderName + "." + shaderType);

        if (!shaderString.empty())
        {
            GLuint shader = glCreateShader(glShaderTypes[i]);
            int len = shaderString.length();
            const char *shaderCStr = shaderString.c_str();
            glShaderSource(shader, 1, (const GLchar **)&shaderCStr, &len);

            glCompileShader(shader);

            GLint compiled;
            glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
            if (compiled == GL_FALSE)
            {
                cout << "Error compiling " << shaderType << " shader." << endl;
                int infoLogLen = 0;
                int charsWritten = 0;
                GLchar *infoLog;

                glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

                if (infoLogLen > 0)
                {
                    infoLog = new GLchar[infoLogLen];
                    glGetShaderInfoLog(shader, infoLogLen, &charsWritten, infoLog);
                    cout << "InfoLog : " << endl
                         << infoLog << endl;
                    delete[] infoLog;
                }

                exit(-1);
            }
            else
            {
                compiledShaders.push_back(shader);
                cout << "Compiled a " << shaderType << " shader" << endl;
            }
        }
    }

    this->shaderProgram = glCreateProgram();

    for (GLuint shaderID : compiledShaders)
    {
        glAttachShader(this->shaderProgram, shaderID);
    }

    glBindAttribLocation(this->shaderProgram, 0, "vertex");
    glBindAttribLocation(this->shaderProgram, 1, "normal");
    glBindAttribLocation(this->shaderProgram, 2, "tex_coords");

    glLinkProgram(this->shaderProgram);
    GLint linked;
    glGetProgramiv(this->shaderProgram, GL_LINK_STATUS, (GLint *)&linked);
    if (linked == GL_FALSE)
    {
        cout << "Failed to link shader." << endl;

        GLint maxLength;
        glGetProgramiv(this->shaderProgram, GL_INFO_LOG_LENGTH, &maxLength);
        if (maxLength > 0)
        {
            char *pLinkInfoLog = new char[maxLength];
            glGetProgramInfoLog(this->shaderProgram, maxLength, &maxLength, pLinkInfoLog);
            cout << pLinkInfoLog << endl;
            delete[] pLinkInfoLog;
        }

        for (GLuint shaderID : compiledShaders)
        {
            glDetachShader(this->shaderProgram, shaderID);
            glDeleteShader(shaderID);
        }

        glDeleteProgram(this->shaderProgram);
        this->shaderProgram = 0;

        exit(-1);
    }
    else {
        cout << "Shader program " << shaderName << " compiled and linked." << endl;
    }
}

void Shader::BindShader() {
    glUseProgram(this->shaderProgram);
}

void Shader::SetupUniforms() {

}