//
//  main.cpp
//  3D Object Drawing
//
//  Created by Nazirul Hasan on 4/9/23.
//

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader.h"
#include "camera.h"
#include "basic_camera.h"
#include "pointLight.h"
#include "cube.h"
#include "stb_image.h"
#include "sphere.h"
#include "Sphere2.h"
#include "Pyramid.h"

#include <iostream>

using namespace std;

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
//void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax);

//void bed(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether);
void bed(Cube &cube, Shader& lightingShaderWithTexture, Shader & lightingShader, glm::mat4 alTogether);
void room(Cube &cube, Shader & lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
//void rightwall(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 alTogether);
void table(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
void Fan(Cube& cube, Shader& lightingShader, glm::mat4 alTogether);
void chair(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
void stair(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
void pond(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether);
void KhanJahanAliHall(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture);
void board(Cube &cube,Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture);
//void drawFish(Shader& lightingShader, glm::mat4 model);
//void road(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture, unsigned int& cVAO);
void road(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture);
// settings

const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

float moveFish = 5, sign = 1;

// modelling transform
float rotateAngle_X = 0.0;
float rotateAngle_Y = 0.0;
float rotateAngle_Z = 0.0;
float rotateAxis_X = 0.0;
float rotateAxis_Y = 0.0;
float rotateAxis_Z = 1.0;
float translate_X = 0.0;
float translate_Y = 0.0;
float translate_Z = 0.0;
float scale_X = 1.0;
float scale_Y = 1.0;
float scale_Z = 1.0;

float eyeX = 8.0f, eyeY = 7.0f, eyeZ = 90.0f;
float lookAtX = 0.0, lookAtY = 0.0, lookAtZ = 0.0;
// camera
Camera camera(glm::vec3(eyeX, eyeY, eyeZ));
float lastX = SCR_WIDTH / 2.0f;
float lastY = SCR_HEIGHT / 2.0f;
bool firstMouse = true;

glm::vec3 V = glm::vec3(0.0f, 1.0f, 0.0f);
BasicCamera basic_camera(eyeX, eyeY, eyeZ, lookAtX, lookAtY, lookAtZ, V);
float rotateFan = 0;
Sphere2 *sphere;

// positions of the point lights
glm::vec3 pointLightPositions[] = {
    glm::vec3(3.50f,  5.50f,  20.0f), /// Spot
    glm::vec3(-1.0f,  1.5f,  0.0f),   /// Point
    glm::vec3(3.0f,  30,  0.5f), /// Sun
    //glm::vec3(-1.5f,  -1.5f,  0.0f)
};
PointLight pointlight1(

    pointLightPositions[0].x, pointLightPositions[0].y, pointLightPositions[0].z,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.8f, 0.8f, 0.8f,     // ambient
    0.8f, 0.8f, 0.8f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    1       // light number
);
PointLight pointlight2(

    pointLightPositions[1].x, pointLightPositions[1].y, pointLightPositions[1].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    1.0f, 1.0f, 1.0f,       // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    2       // light number
);

PointLight pointlight3(

    pointLightPositions[2].x, pointLightPositions[2].y, pointLightPositions[2].z,  // position
    0.7f, 0.7f, 0.7f,     // ambient
    0.7f, 0.7f, 0.7f,      // diffuse
    0.7f, 0.7f, 0.7f,   // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    100       // light number
);


PointLight pointlight4(

    12.5-0.08-0.5, 2.8+2, 54-0.02-0.5-25,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.8f, 0.8f, 0.8f,     // ambient
    0.8f, 0.8f, 0.8f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    4       // light number
);



PointLight pointlight5(

    12.5 - 0.08 - 0.5, 2.8 + 2, 54 - 0.02 - 0.5 - 15,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.8f, 0.8f, 0.8f,     // ambient
    0.8f, 0.8f, 0.8f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    5       // light number
);

PointLight pointlight6(

    12.5 - 0.08 - 0.5, 2.8 + 2, 54 - 0.02 - 0.5 - 5,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.8f, 0.8f, 0.8f,     // ambient
    0.8f, 0.8f, 0.8f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    6       // light number
);


PointLight pointlight7(

    12.5 - 0.08 - 0.5, 2.8 + 2, 54 - 0.02 - 0.5 +5,  // position
    //1.0f, 1.0f, 1.0f,     // ambient
    //1.0f, 1.0f, 1.0f,      // diffuse
    //1.0f, 1.0f, 1.0f,        // specular
    0.8f, 0.8f, 0.8f,     // ambient
    0.8f, 0.8f, 0.8f,      // diffuse
    1.0f, 1.0f, 1.0f,        // specular
    1.0f,   //k_c
    0.09f,  //k_l
    0.032f, //k_q
    7       // light number
);

// light settings
bool onOffPointToggle = true;
bool onOffSpotToggle = true;
bool onOffDirectToggle = true;
bool ambientToggle = true;
bool diffuseToggle = true;
bool specularToggle = true;

float doorangle = 0.0f;

// timing
float deltaTime = 0.0f;    // time between current frame and last frame
float lastFrame = 0.0f;


    unsigned int cubeVAO, cubeVBO, cubeEBO;
    unsigned int cVAO, cVBO, cEBO;

    void drawCube(unsigned int& cubeVAO, Shader& lightingShader, glm::mat4 model = glm::mat4(1.0f), float r = 1.0f, float g = 1.0f, float b = 1.0f)
    {
        lightingShader.use();

        lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
        lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
        lightingShader.setVec3("material.specular", glm::vec3(0.5f, 0.5f, 0.5f));
        lightingShader.setFloat("material.shininess", 32.0f);

        lightingShader.setMat4("model", model);

        glBindVertexArray(cubeVAO);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    }
    void drawCylinder(unsigned int& cVAO, Shader& lightingShader, glm::mat4 model, float r, float g, float b)
    {
        lightingShader.use();

        lightingShader.setVec3("material.ambient", glm::vec3(r, g, b));
        lightingShader.setVec3("material.diffuse", glm::vec3(r, g, b));
        lightingShader.setVec3("material.specular", glm::vec3(1.0f, 1.0f, 1.0f));
        lightingShader.setFloat("material.shininess", 60.0f);

        lightingShader.setMat4("model", model);

        glBindVertexArray(cVAO);
        glDrawElements(GL_TRIANGLES, 120, GL_UNSIGNED_INT, 0);
    }

class Curve
{
public:
    vector<float> cntrlPoints;
    vector <float> coordinates;
    vector <float> normals;
    vector <int> indices;
    vector <float> vertices;
    const double pi = 3.14159265389;
    const int nt = 40;
    const int ntheta = 20;
    Curve(vector<float> &tmp)
    {
        this->cntrlPoints = tmp;
        this->fishVAO = hollowBezier(cntrlPoints.data(), ((unsigned int)cntrlPoints.size() / 3) - 1);
        cout << cntrlPoints.size() << endl;
        cout << coordinates.size() << endl;
        cout << normals.size() << endl;
        cout << indices.size() << endl;
        cout << vertices.size() << endl;
    }
    ~Curve()
    {
        glDeleteVertexArrays(1, &fishVAO);
        glDeleteVertexArrays(1, &bezierVAO);
        glDeleteBuffers(1, &bezierVBO);
        glDeleteBuffers(1, &bezierEBO);
    }
    void draw(Shader& lightingShader, glm::mat4 model)
    {
        /// Fish
        lightingShader.use();
        lightingShader.setMat4("model", model);
        lightingShader.setVec3("material.ambient", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.diffuse", glm::vec3(1.0f, 0.6f, 0.0f));
        lightingShader.setVec3("material.specular", glm::vec3(1.0f, 1.0f, 1.0f));
        lightingShader.setFloat("material.shininess", 32.0f);

        glBindVertexArray(fishVAO);
        glDrawElements(GL_TRIANGLES,                    // primitive type
            (unsigned int)indices.size(),          // # of indices
            GL_UNSIGNED_INT,                 // data type
            (void*)0);                       // offset to indices

        // unbind VAO
        glBindVertexArray(0);
        /// End Fish
    }
private:
    unsigned int fishVAO;
    unsigned int bezierVAO;
    unsigned int bezierVBO;
    unsigned int bezierEBO;


    unsigned int drawControlPoints()
    {
        unsigned int controlPointVAO;
        unsigned int controlPointVBO;

        glGenVertexArrays(1, &controlPointVAO);
        glGenBuffers(1, &controlPointVBO);

        glBindVertexArray(controlPointVAO);

        glBindBuffer(GL_ARRAY_BUFFER, controlPointVBO);
        glBufferData(GL_ARRAY_BUFFER, (unsigned int)cntrlPoints.size() * sizeof(float), cntrlPoints.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        return controlPointVAO;
    }

    long long nCr(int n, int r)
    {
        if (r > n / 2)
            r = n - r; // because C(n, r) == C(n, n - r)
        long long ans = 1;
        int i;

        for (i = 1; i <= r; i++)
        {
            ans *= n - r + i;
            ans /= i;
        }

        return ans;
    }
    void BezierCurve(double t, float xy[2], GLfloat ctrlpoints[], int L)
    {
        double y = 0;
        double x = 0;
        t = t > 1.0 ? 1.0 : t;
        for (int i = 0; i < L + 1; i++)
        {
            long long ncr = nCr(L, i);
            double oneMinusTpow = pow(1 - t, double(L - i));
            double tPow = pow(t, double(i));
            double coef = oneMinusTpow * tPow * ncr;
            x += coef * ctrlpoints[i * 3];
            y += coef * ctrlpoints[(i * 3) + 1];

        }
        xy[0] = float(x);
        xy[1] = float(y);
    }
    unsigned int hollowBezier(GLfloat ctrlpoints[], int L)
    {
        int i, j;
        float x, y, z, r;                //current coordinates
        float theta;
        float nx, ny, nz, lengthInv;    // vertex normal


        const float dtheta = 2 * pi / ntheta;        //angular step size

        float t = 0;
        float dt = 1.0 / nt;
        float xy[2];

        for (i = 0; i <= nt; ++i)              //step through y
        {
            BezierCurve(t, xy, ctrlpoints, L);
            r = xy[0];
            y = xy[1];
            theta = 0;
            t += dt;
            lengthInv = 1.0 / r;

            for (j = 0; j <= ntheta; ++j)
            {
                double cosa = cos(theta);
                double sina = sin(theta);
                z = r * cosa;
                x = r * sina;

                coordinates.push_back(x);
                coordinates.push_back(y);
                coordinates.push_back(z);

                // normalized vertex normal (nx, ny, nz)
                // center point of the circle (0,y,0)
                nx = (x - 0) * lengthInv;
                ny = (y - y) * lengthInv;
                nz = (z - 0) * lengthInv;

                normals.push_back(nx);
                normals.push_back(ny);
                normals.push_back(nz);

                theta += dtheta;
            }
        }
        // generate index list of triangles
        // k1--k1+1
        // |  / |
        // | /  |
        // k2--k2+1

        int k1, k2;
        for (int i = 0; i < nt; ++i)
        {
            k1 = i * (ntheta + 1);     // beginning of current stack
            k2 = k1 + ntheta + 1;      // beginning of next stack

            for (int j = 0; j < ntheta; ++j, ++k1, ++k2)
            {
                // k1 => k2 => k1+1
                indices.push_back(k1);
                indices.push_back(k2);
                indices.push_back(k1 + 1);

                // k1+1 => k2 => k2+1
                indices.push_back(k1 + 1);
                indices.push_back(k2);
                indices.push_back(k2 + 1);
            }
        }

        size_t count = coordinates.size();
        for (int i = 0; i < count; i += 3)
        {
            //cout << count << ' ' << i + 2 << endl;
            vertices.push_back(coordinates[i]);
            vertices.push_back(coordinates[i + 1]);
            vertices.push_back(coordinates[i + 2]);

            vertices.push_back(normals[i]);
            vertices.push_back(normals[i + 1]);
            vertices.push_back(normals[i + 2]);
        }

        glGenVertexArrays(1, &bezierVAO);
        glBindVertexArray(bezierVAO);

        // create VBO to copy vertex data to VBO
        glGenBuffers(1, &bezierVBO);
        glBindBuffer(GL_ARRAY_BUFFER, bezierVBO);           // for vertex data
        glBufferData(GL_ARRAY_BUFFER,                   // target
            (unsigned int)vertices.size() * sizeof(float), // data size, # of bytes
            vertices.data(),   // ptr to vertex data
            GL_STATIC_DRAW);                   // usage

        // create EBO to copy index data
        glGenBuffers(1, &bezierEBO);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bezierEBO);   // for index data
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,           // target
            (unsigned int)indices.size() * sizeof(unsigned int),             // data size, # of bytes
            indices.data(),               // ptr to index data
            GL_STATIC_DRAW);                   // usage

        // activate attrib arrays
        glEnableVertexAttribArray(0);
        glEnableVertexAttribArray(1);

        // set attrib arrays with stride and offset
        int stride = 24;     // should be 24 bytes
        glVertexAttribPointer(0, 3, GL_FLOAT, false, stride, (void*)0);
        glVertexAttribPointer(1, 3, GL_FLOAT, false, stride, (void*)(sizeof(float) * 3));

        // unbind VAO, VBO and EBO
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

        return bezierVAO;
    }

};

vector<float>Fish = {
-0.0100, 1.9950, 5.1000,
-0.0550, 1.9800, 5.1000,
-0.0950, 1.9350, 5.1000,
-0.1500, 1.8250, 5.1000,
-0.2250, 1.5900, 5.1000,
-0.2550, 1.3450, 5.1000,
-0.2050, 1.1700, 5.1000,
-0.1400, 1.0050, 5.1000,
-0.0400, 0.8600, 5.1000,
0.0400, 0.7300, 5.1000,
0.1300, 0.6350, 5.1000,
0.2400, 0.5050, 5.1000,
};
Cube* tmp, * roomwindow, * roomfloor, * grass, *roomdoor, *walltex, *pondtex, *road1, *divider, *khaja;
Curve* fis;
int ind = 0;

Pyramid *pyramid;

void useShaderProgram(Shader& lightingShaderWithTexture)
{
    lightingShaderWithTexture.use();
    pointlight1.setUpPointLight(lightingShaderWithTexture);
    // point light 2
    pointlight2.setUpPointLight(lightingShaderWithTexture);
    // point light 3
    pointlight3.setUpPointLight(lightingShaderWithTexture);
    pointlight4.setUpPointLight(lightingShaderWithTexture);

    pointlight5.setUpPointLight(lightingShaderWithTexture);

    pointlight6.setUpPointLight(lightingShaderWithTexture);
    pointlight7.setUpPointLight(lightingShaderWithTexture);
}
int main()
{

    // glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "CSE 4208: Computer Graphics Laboratory", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    //glfwSetKeyCallback(window, key_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    // tell GLFW to capture our mouse
    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_HIDDEN);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }


    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // build and compile our shader zprogram
    // ------------------------------------
    
    const int parts = 50;
    const float pi = 3.1415926535;
    const float angle = pi * 2.0f / parts;
    float points[200000]{}, radius = 1.0f;

    points[ind++] = 0.0f;
    points[ind++] = 0.0f;
    points[ind++] = 0.0f;
    for (int i = 1; i <= parts; i++) {
        points[ind++] = radius * sin(angle * i);
        points[ind++] = -radius * cos(angle * i);
        points[ind++] = 0.0f;
    }

    for (float r = radius - 0.005f, z = 0.005f; r > 0.0f; r -= 0.005f, z += 0.005f)
    {
        for (int i = 1; i <= parts + 1; i++) {
            points[ind++] = (r + 0.005) * sin(angle * i);
            points[ind++] = -(r + 0.005) * cos(angle * i);
            points[ind++] = z - 0.005f;

            points[ind++] = r * sin(angle * i);
            points[ind++] = -r * cos(angle * i);
            points[ind++] = z;
        }
    }
    for (float r = radius - 0.005f, z = -0.005f; r > 0.0f; r -= 0.005f, z -= 0.005f)
    {
        for (int i = 1; i <= parts + 1; i++) {
            points[ind++] = (r + 0.005) * sin(angle * i);
            points[ind++] = -(r + 0.005) * cos(angle * i);
            points[ind++] = z + 0.005f;

            points[ind++] = r * sin(angle * i);
            points[ind++] = -r * cos(angle * i);
            points[ind++] = z;
        }
    }

    /// Sphere
    unsigned int VBOCL, shpareVAO;
    glGenVertexArrays(1, &shpareVAO);
    glGenBuffers(1, &VBOCL);
    glBindVertexArray(shpareVAO);
    glBindBuffer(GL_ARRAY_BUFFER, VBOCL);
    glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);
    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    //cyllinder

    float ver_arr[] = {

    1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f,
    0.809017f, 1.0f, 0.587785f, 0.809017f, 1.0f, 0.587785f,
    0.309017f, 1.0f, 0.951057f, 0.309017f, 1.0f, 0.951057f,
    -0.309017f, 1.0f, 0.951057f, -0.309017f, 1.0f, 0.951057f,
    -0.809017f, 1.0f, 0.587785f, -0.809017f, 1.0f, 0.587785f,
    -1.0f, 1.0f, 1.22465e-16f, -1.0f, 1.0f, 1.22465e-16f,
    -0.809017f, 1.0f, -0.587785f, -0.809017f, 1.0f, -0.587785f,
    -0.309017f, 1.0f, -0.951057f, -0.309017f, 1.0f, -0.951057f,
    0.309017f, 1.0f, -0.951057f, 0.309017f, 1.0f, -0.951057f,
    0.809017f, 1.0f, -0.587785f, 0.809017f, 1.0f, -0.587785f,

    1.0f, -1.0f, 0.0f, 1.0f, -1.0f, 0.0f,
    0.809017f, -1.0f, 0.587785f, 0.809017f, -1.0f, 0.587785f,
    0.309017f, -1.0f, 0.951057f, 0.309017f, -1.0f, 0.951057f,
    -0.309017f, -1.0f, 0.951057f, -0.309017f, -1.0f, 0.951057f,
    -0.809017f, -1.0f, 0.587785f, -0.809017f, -1.0f, 0.587785f,
    -1.0f, -1.0f, 1.22465e-16f, -1.0f, -1.0f, 1.22465e-16f,
    -0.809017f, -1.0f, -0.587785f, -0.809017f, -1.0f, -0.587785f,
    -0.309017f, -1.0f, -0.951057f, -0.309017f, -1.0f, -0.951057f,
    0.309017f, -1.0f, -0.951057f, 0.309017f, -1.0f, -0.951057f,
    0.809017f, -1.0f, -0.587785f, 0.809017f, -1.0f, -0.587785f,


    1.0f, -1.0f, 0.0f, 1.0f, -1.0f, 0.0f,
    0.809017f, -1.0f, 0.587785f, 0.809017f, -1.0f, 0.587785f,
    0.309017f, -1.0f, 0.951057f, 0.309017f, -1.0f, 0.951057f,
    -0.309017f, -1.0f, 0.951057f, -0.309017f, -1.0f, 0.951057f,
    -0.809017f, -1.0f, 0.587785f, -0.809017f, -1.0f, 0.587785f,
    -1.0f, -1.0f, 1.22465e-16f, -1.0f, -1.0f, 1.22465e-16f,
    -0.809017f, -1.0f, -0.587785f, -0.809017f, -1.0f, -0.587785f,
    -0.309017f, -1.0f, -0.951057f, -0.309017f, -1.0f, -0.951057f,
    0.309017f, -1.0f, -0.951057f, 0.309017f, -1.0f, -0.951057f,
    0.809017f, -1.0f, -0.587785f, 0.809017f, -1.0f, -0.587785f,

    1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f,
    0.809017f, 1.0f, 0.587785f, 0.809017f, 1.0f, 0.587785f,
    0.309017f, 1.0f, 0.951057f, 0.309017f, 1.0f, 0.951057f,
    -0.309017f, 1.0f, 0.951057f, -0.309017f, 1.0f, 0.951057f,
    -0.809017f, 1.0f, 0.587785f, -0.809017f, 1.0f, 0.587785f,
    -1.0f, 1.0f, 1.22465e-16f, -1.0f, 1.0f, 1.22465e-16f,
    -0.809017f, 1.0f, -0.587785f, -0.809017f, 1.0f, -0.587785f,
    -0.309017f, 1.0f, -0.951057f, -0.309017f, 1.0f, -0.951057f,
    0.309017f, 1.0f, -0.951057f, 0.309017f, 1.0f, -0.951057f,
    0.809017f, 1.0f, -0.587785f, 0.809017f, 1.0f, -0.587785f,


    0.0,-1.0,0.0, 0.0,-1.0,0.0,
    0.0,1.0,0.0, 0.0,1.0,0.0

    };

    unsigned int ind_arr[] = {
        0, 11, 1,
        11, 0, 10,
        1, 12, 2,
        12, 1, 11,
        2, 13, 3,
        13, 2, 12,
        3, 14, 4,
        14, 3, 13,
        4, 15, 5,
        15, 4, 14,
        5, 16, 6,
        16, 5, 15,
        6, 17, 7,
        17, 6, 16,
        7, 18, 8,
        18, 7, 17,
        8, 19, 9,
        19, 8, 18,
        9, 10, 0,
        10, 9, 19,

        40,20,21,
        40,21,22,
        40,22,23,
        40,23,24,
        40,24,25,
        40,25,26,
        40,26,27,
        40,27,28,
        40,28,29,
        40,29,20,

        41,30,31,
        41,31,32,
        41,32,33,
        41,33,34,
        41,34,35,
        41,35,36,
        41,36,37,
        41,37,38,
        41,38,39,
        41,39,30


        /*
        21,10,11,
        21,11,12,
        21,12,13,
        21,13,14,
        21,14,15,
        21,15,16,
        21,16,17,
        21,17,18,
        21,18,19,
        21,19,10*/

    };

    glGenVertexArrays(1, &cVAO);
    glGenBuffers(1, &cVBO);
    glGenBuffers(1, &cEBO);

    glBindVertexArray(cVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(ver_arr), ver_arr, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(ind_arr), ind_arr, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);
    //end of cylingder

    float cube_vertices[] = {
        // positions      // normals
        0.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 0.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 0.0f, -1.0f,

        1.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 0.0f,
        1.0f, 1.0f, 1.0f, 1.0f, 0.0f, 0.0f,

        0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 0.0f, 1.0f,

        0.0f, 0.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 1.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 1.0f, 0.0f, -1.0f, 0.0f, 0.0f,
        0.0f, 0.0f, 0.0f, -1.0f, 0.0f, 0.0f,

        1.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,
        1.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 0.0f, 0.0f, 1.0f, 0.0f,
        0.0f, 1.0f, 1.0f, 0.0f, 1.0f, 0.0f,

        0.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 0.0f, 0.0f, -1.0f, 0.0f,
        1.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f,
        0.0f, 0.0f, 1.0f, 0.0f, -1.0f, 0.0f
    };
    unsigned int cube_indices[] = {
        0, 3, 2,
        2, 1, 0,

        4, 5, 7,
        7, 6, 4,

        8, 9, 10,
        10, 11, 8,

        12, 13, 14,
        14, 15, 12,

        16, 17, 18,
        18, 19, 16,

        20, 21, 22,
        22, 23, 20
    };

    glGenVertexArrays(1, &cubeVAO);
    glGenBuffers(1, &cubeVBO);
    glGenBuffers(1, &cubeEBO);

    glBindVertexArray(cubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(cube_vertices), cube_vertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cube_indices), cube_indices, GL_STATIC_DRAW);


    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // vertex normal attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)12);
    glEnableVertexAttribArray(1);

    // second, configure the light's VAO (VBO stays the same; the vertices are the same for the light object which is also a 3D cube)
    unsigned int lightCubeVAO;
    glGenVertexArrays(1, &lightCubeVAO);
    glBindVertexArray(lightCubeVAO);

    glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cubeEBO);
    // note that we update the lamp's position attribute's stride to reflect the updated buffer data
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);


    Shader lightingShaderWithTexture("vertexShaderForPhongShadingWithTexture.vs", "fragmentShaderForPhongShadingWithTexture.fs");
    Shader lightingShader("vertexShaderForPhongShading.vs", "fragmentShaderForPhongShading.fs");
    Shader ourShader("vertexShader.vs", "fragmentShader.fs");

    //Sphere sphere = Sphere();

    string diffuseMapPath = "";
    string specularMapPath = "";

    diffuseMapPath = "Wall2.png";
    specularMapPath = "container2_specular.png";


    unsigned int diffMap = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube = Cube(diffMap, specMap, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);


    diffuseMapPath = "Tile2.jpg";
    specularMapPath = "whiteBackground.png";


    unsigned int diffMap2 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap2 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube2 = Cube(diffMap2, specMap2, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    tmp = &cube2;

    diffuseMapPath = "Window.png";
    specularMapPath = "WindowSpec.jpg";

    unsigned int diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    unsigned int specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube3 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    roomwindow = &cube3;

    diffuseMapPath = "roomFloor.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube4 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    roomfloor = &cube4;

    diffuseMapPath = "grass.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube5 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    grass = &cube5;

    diffuseMapPath = "door.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube6 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    roomdoor = &cube6;
    
    diffuseMapPath = "Walltex.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube7 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    walltex = &cube7;

    diffuseMapPath = "pondtex.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube8 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    pondtex = &cube8;


    /*diffuseMapPath = "pondtex.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube9 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    pondtex = &cube9;*/


    //Sphere x = Sphere();
    //sphere = &x;
    //
    diffuseMapPath = "lamp.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Sphere2 sp(1.0, 36, 18, glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(1.0f, 0.0f, 0.0f), glm::vec3(0.5f, 0.5f, 0.5f), 32.0f, diffMap3, specMap3, 0, 1, 0, 1);
    sp.setDefaults();
    sp.setTexture(diffMap3,specMap3);
    sphere = &sp;

    Sphere sun(1.0f,36,18,glm::vec3(1.0f, 1.0f, 1.0f),glm::vec3(1.0f,1.0f,1.0f),glm::vec3(1.0f, 1.0f,1.0f),32.0f);

    diffuseMapPath = "road.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube10 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    road1 = &cube10;

    diffuseMapPath = "divider.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube11 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    divider = &cube11;

    diffuseMapPath = "khaja.png";
    specularMapPath = "WindowSpec.jpg";

    diffMap3 = loadTexture(diffuseMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    specMap3 = loadTexture(specularMapPath.c_str(), GL_REPEAT, GL_REPEAT, GL_LINEAR_MIPMAP_LINEAR, GL_LINEAR);
    Cube cube12 = Cube(diffMap3, specMap3, 32.0f, 0.0f, 0.0f, 1.0f, 1.0f);
    khaja = &cube12;
     
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    Pyramid abc("tree.png");
    pyramid = &abc;

    Curve fish(Fish);
    fis = &fish;

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // per-frame time logic
        // --------------------
        float currentFrame = static_cast<float>(glfwGetTime());
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        // input
        // -----
        processInput(window);

        // render
        // ------
        if (pointlight3.isOn())
        {
            glClearColor(0.5294f, 0.8078f, 0.9804f, 0.08f);
        //    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
        }
        else
        {
            glClearColor(0.4f, 0.4f, 0.2f, 1.0f);
          //  ourShader.setVec3("color", glm::vec3(0.2f, 0.2f, 0.2f));
        }
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // be sure to activate shader when setting uniforms/drawing objects
        lightingShaderWithTexture.use();
        lightingShaderWithTexture.setVec3("viewPos", camera.Position);

        // pass projection matrix to shader (note that in this case it could change every frame)
        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
        //glm::mat4 projection = glm::ortho(-2.0f, +2.0f, -1.5f, +1.5f, 0.1f, 100.0f);
        lightingShaderWithTexture.setMat4("projection", projection);

        // camera/view transformation
        glm::mat4 view = camera.GetViewMatrix();
        //glm::mat4 view = basic_camera.createViewMatrix();
        lightingShaderWithTexture.setMat4("view", view);

        lightingShader.use();
        lightingShader.setVec3("viewPos", camera.Position);
        lightingShader.setMat4("projection", projection);
        lightingShader.setMat4("view", view);
        

        //pointlight1.setUpPointLight(lightingShader);
        //pointlight2.setUpPointLight(lightingShader);
        //pointlight3.setUpPointLight(lightingShader);

        // Modelling Transformation
        glm::mat4 identityMatrix = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
        glm::mat4 translateMatrix, rotateXMatrix, rotateYMatrix, rotateZMatrix, scaleMatrix, model;
        translateMatrix = glm::translate(identityMatrix, glm::vec3(translate_X, translate_Y, translate_Z));
        rotateXMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_X), glm::vec3(1.0f, 0.0f, 0.0f));
        rotateYMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Y), glm::vec3(0.0f, 1.0f, 0.0f));
        rotateZMatrix = glm::rotate(identityMatrix, glm::radians(rotateAngle_Z), glm::vec3(0.0f, 0.0f, 1.0f));
        scaleMatrix = glm::scale(identityMatrix, glm::vec3(scale_X, scale_Y, scale_Z));
        model = translateMatrix * rotateXMatrix * rotateYMatrix * rotateZMatrix * scaleMatrix;
 
        //glm::mat4 modelMatrixForContainer = glm::mat4(1.0f);
        //modelMatrixForContainer = glm::translate(model, glm::vec3(-0.45f, -0.4f, -2.8f));
        //cube.drawCubeWithTexture(lightingShaderWithTexture, modelMatrixForContainer);

        useShaderProgram(lightingShader);
        glm::mat4 tmp = model, tempt, save;
        tempt = glm::translate(identityMatrix, glm::vec3(-6, -1, -7.5));
        tmp = model * tempt;
        pond(cube, lightingShaderWithTexture, lightingShader, tmp* glm::translate(identityMatrix, glm::vec3(-35, 0, 0)));
        save = tmp;
        //Fan(cube, lightingShader, tmp);
        float qqq = 6.0;
        float ppp = 0;
        for (int i = 0; i < 3; i++)
        {
            room(cube, lightingShaderWithTexture, lightingShader, tmp);
            ppp = 0;
            for (int i = 0; i < 6; i++)
            {
                tmp *= glm::translate(identityMatrix, glm::vec3(6, 0, 0));
                ppp += 6.0f;
                room(cube, lightingShaderWithTexture, lightingShader, tmp);
            }
            tmp = save;
            tmp *= glm::translate(identityMatrix, glm::vec3(0, qqq, 0));
            qqq += 6.0;
        }

        glm::mat4 rYM = glm::rotate(identityMatrix, glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        tempt = glm::translate(identityMatrix, glm::vec3(ppp+10, -1, 2));
        tmp = model * tempt * rYM;
        save = tmp;
        qqq = 6.0;
        for (int i = 0; i < 3; i++)
        {
            room(cube, lightingShaderWithTexture, lightingShader, tmp);
            for (int i = 0; i < 6; i++)
            {
                tmp *= glm::translate(identityMatrix, glm::vec3(6, 0, 0));
                room(cube, lightingShaderWithTexture, lightingShader, tmp);
            }
            tmp = save;
            tmp *= glm::translate(identityMatrix, glm::vec3(0, qqq, 0));
            qqq += 6.0;
        }
        //stair(cube, lightingShaderWithTexture, lightingShader, save);
        //bed(cube, lightingShaderWithTexture, lightingShader, modelMatrixForContainer);

        useShaderProgram(lightingShaderWithTexture);
        

        board(cube, lightingShader, model * glm::translate(glm::mat4(1), glm::vec3(4, 0, 26)), lightingShaderWithTexture);
        
        KhanJahanAliHall(lightingShader, model*glm::translate(glm::mat4(1),glm::vec3(12.5, 2, 54)), lightingShaderWithTexture);
        road(lightingShader, model * glm::translate(glm::mat4(1), glm::vec3(12, -0.95, 50)), lightingShaderWithTexture);


        lightingShaderWithTexture.use();
        glm::mat4 floormat = glm::translate(identityMatrix, glm::vec3(-100, -1.2, -100));
        glm::mat4 floorscale = glm::scale(identityMatrix, glm::vec3(1000, 0.2, 1000));
        glm::mat4 floormodel = model * floormat * floorscale;
        grass->drawCubeWithTexture(lightingShaderWithTexture, floormodel);

        // also draw the lamp object(s)
        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        // we now draw as many light bulbs as we have point lights.
        //for (unsigned int i = 0; i < 4; i++)
        //{
        //    model = glm::mat4(1.0f);
        //    model = glm::translate(model, pointLightPositions[i]);
        //    model = glm::scale(model, glm::vec3(0.2f)); // Make it a smaller cube
        //    cube.drawCube(ourShader, model, 0.8f, 0.8f, 0.8f);
        //}

        glBindVertexArray(shpareVAO);
        for (unsigned int i = 0; i < 3; i++)
        {
            model = glm::mat4(1.0f);
            model = glm::translate(model, pointLightPositions[i]);
            model = glm::scale(model, glm::vec3(0.2f)); // Make it a smaller cube
            tempt = glm::translate(identityMatrix, glm::vec3(-20, -1, -30));
            ourShader.setMat4("model", model * tempt);
            if (i == 0)
            {
                if (pointlight1.isOn())
                    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
                else
                    ourShader.setVec3("color", glm::vec3(0.2f, 0.2f, 0.2f));
            }
            if (i == 1)
            {
                if (pointlight2.isOn())
                    ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
                else
                    ourShader.setVec3("color", glm::vec3(0.2f, 0.2f, 0.2f));
            }
            if (i == 2)
            {
                //model = glm::mat4(1.0f);
                //model = glm::translate(model, pointLightPositions[i]);
                //model = glm::scale(model, glm::vec3(2.0f)); // Make it a smaller cube
                //tempt = glm::translate(identityMatrix, glm::vec3(-20, -1, -30));

                //    useShaderProgram(lightingShaderWithTexture);
                //    lightingShaderWithTexture.use();
                //    sp.drawSphereWithTexture(lightingShaderWithTexture, model);
                //    useShaderProgram(lightingShader);
                if (pointlight3.isOn())
                {
                    sun.setColor(glm::vec3(1.0f, 0.8f, 0.0f));
                    sun.drawSphere(lightingShader, model * tempt * glm::scale(glm::mat4(1.0f), glm::vec3(6, 6, 6)));
                    //sphere->drawSphere(lightingShader, model * tempt);
                    //ourShader.setVec3("color", glm::vec3(0.8f, 0.8f, 0.8f));
                }
                else
                {
                    //sun.setColor(glm::vec3(0.1, 0.1, 0.1));
                    //sun.drawSphere(lightingShader, model* tempt*glm::scale(glm::mat4(1.0f), glm::vec3(6, 6, 6)));
                    glClearColor(0.6f, 0.4f, 0.2f, 1.0f);
                    //ourShader.setVec3("color", glm::vec3(0.2f, 0.2f, 0.2f));
                }
            }

            //glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, ind);
        }

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------


    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}
void board(Cube &cube, Shader &lightingShader, glm::mat4 alTogether, Shader &lightingShaderWithTexture)
{
    float height = 3;
    float width = 15;
    float pur = 0.3;

    float pilheight = 10;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(0.0, 0.0 + pilheight, 0.0));
    scale = glm::scale(model, glm::vec3(width, height, pur));
    glm::mat4 rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rot * scale;
    khaja->drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.3, 0.6,0.9);


    // piller 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.5, pilheight, 1));
    rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rot * scale;
    drawCylinder(cVAO, lightingShader, model, 0.3, 0.6, 0.9);

    float treeheight = 5;

    /// Tree

    for (float dist = 8.0; dist <= 35; dist += 10)
    {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate = glm::translate(model, glm::vec3(-3, 0.0, dist));
        scale = glm::scale(model, glm::vec3(0.5, treeheight, 1));
        rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        model = alTogether * translate * rot * scale;
        drawCylinder(cVAO, lightingShader, model, 0.5, 0.5, 0.5);
        useShaderProgram(lightingShaderWithTexture);
        for (float i = 0.0; i <= 5; i += 1.2)
        {
            model = glm::mat4(1.0f);
            translate = glm::mat4(1.0f);
            translate2 = glm::mat4(1.0f);
            scale = glm::mat4(1.0f);
            translate = glm::translate(model, glm::vec3(-3.0, treeheight-i, dist));
            scale = glm::scale(model, glm::vec3(5, 3.2, 5));
            model = alTogether * translate * scale;
            pyramid->draw(lightingShaderWithTexture, model);
        }
    }
    // row
    for (float dist = -3; dist >= -40; dist -= 10)
    {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate = glm::translate(model, glm::vec3(0+dist, 0.0, 38));
        scale = glm::scale(model, glm::vec3(0.5, treeheight, 1));
        rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        model = alTogether * translate * rot * scale;
        drawCylinder(cVAO, lightingShader, model, 0.5, 0.5, 0.5);
        useShaderProgram(lightingShaderWithTexture);
        for (float i = 0.0; i <= 5; i += 1.2)
        {
            model = glm::mat4(1.0f);
            translate = glm::mat4(1.0f);
            translate2 = glm::mat4(1.0f);
            scale = glm::mat4(1.0f);
            translate = glm::translate(model, glm::vec3(0+dist, treeheight - i, 38));
            scale = glm::scale(model, glm::vec3(5, 3.2, 5));
            model = alTogether * translate * scale;
            pyramid->draw(lightingShaderWithTexture, model);
        }
    }
    

    // piller 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(width, 0.0, 0.0));
    scale = glm::scale(model, glm::vec3(0.5, pilheight, 1));
    rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rot * scale;
    drawCylinder(cVAO, lightingShader, model, 0.3, 0.6, 0.9);

    // Tree 2
    for (float dist = 8.0; dist <= 35; dist += 10)
    {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate = glm::translate(model, glm::vec3(width + 3, 0.0, dist));
        scale = glm::scale(model, glm::vec3(0.5, treeheight, 1));
        rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        model = alTogether * translate * rot * scale;
        drawCylinder(cVAO, lightingShader, model, 0.5, 0.5, 0.5);
        useShaderProgram(lightingShaderWithTexture);
        for (float i = 0.0; i <= 5; i += 1.2)
        {
            model = glm::mat4(1.0f);
            translate = glm::mat4(1.0f);
            translate2 = glm::mat4(1.0f);
            scale = glm::mat4(1.0f);
            translate = glm::translate(model, glm::vec3(width +3.0, treeheight - i, dist));
            scale = glm::scale(model, glm::vec3(5, 3.2, 5));
            model = alTogether * translate * scale;
            pyramid->draw(lightingShaderWithTexture, model);
        }
    }

    // row

    for (float dist = -3; dist >= -40; dist -= 10)
    {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate = glm::translate(model, glm::vec3(0 - dist + width, 0.0, 38));
        scale = glm::scale(model, glm::vec3(0.5, treeheight, 1));
        rot = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
        model = alTogether * translate * rot * scale;
        drawCylinder(cVAO, lightingShader, model, 0.5, 0.5, 0.5);
        useShaderProgram(lightingShaderWithTexture);
        for (float i = 0.0; i <= 5; i += 1.2)
        {
            model = glm::mat4(1.0f);
            translate = glm::mat4(1.0f);
            translate2 = glm::mat4(1.0f);
            scale = glm::mat4(1.0f);
            translate = glm::translate(model, glm::vec3(0 - dist + width, treeheight - i, 38));
            scale = glm::scale(model, glm::vec3(5, 3.2, 5));
            model = alTogether * translate * scale;
            pyramid->draw(lightingShaderWithTexture, model);
        }
    }

}


void piller(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture)
{
    float width = 0.2;
    float height = 3;


    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, height, width));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * scale * translate;
    useShaderProgram(lightingShader);
    //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
    drawCylinder(cVAO, lightingShader, model, 0 / 255.0, 84 / 255.0, 147 / 255.0);

}
void KhanJahanAliHall(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture)
{
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 translate3 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotate = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 1, 1));
    translate = glm::translate(model, glm::vec3(-0.5, -2, -12));
    model = alTogether * scale * translate;

    //road(lightingShader, model, lightingShaderWithTexture);
    for (int i = -25; i < 15; i += 10) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        translate3 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(0, 0, i));

        //scale = glm::scale(model, glm::vec3(0.1, 3, 0.1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
        model = alTogether * scale * translate * translate2;

        piller(lightingShader, model, lightingShaderWithTexture);

        translate3 = glm::translate(model, glm::vec3(-0.08, 2.8, -0.02));

        scale = glm::scale(model, glm::vec3(0.3, 0.3, 0.3));
        ////translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
        rotate = glm::rotate(rotate, glm::radians(2.0f), glm::vec3(0, 0, 1));
        model = translate3 * rotate * glm::mat4(1.0f);
        sphere->drawSphereWithTexture(lightingShaderWithTexture, model);
    }
}

void road(Shader& lightingShader, glm::mat4 alTogether, Shader& lightingShaderWithTexture)
{
    float road_width = 11;
    float road_length = 40;
    float road_height = 0.3;


    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(road_width, road_height, road_length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
    road1->drawCubeWithTexture(lightingShaderWithTexture, model);
    
    //devider
    for (int i = (-road_length / 2); i < (road_length / 2) - 5; i++) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(-1, road_height, i + 5));
        scale = glm::scale(model, glm::vec3(0.5, 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, -0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
    }

    for (int i = (-road_length / 2); i < (road_length / 2) - 1 - 5; i++) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(1, road_height, i + 5));
        scale = glm::scale(model, glm::vec3(0.5, 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
    }

    for (int i = -road_length / 2; i < road_length / 2; i++) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(road_width / 2, 0, i));
        scale = glm::scale(model, glm::vec3(0.5, road_height + 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0, 0, -0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
    }

    for (int i = -road_length / 2; i < road_length / 2; i++) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(-road_width / 2, 0, i));
        scale = glm::scale(model, glm::vec3(0.5, road_height + 0.5, 1));
        translate = glm::translate(model, glm::vec3(-1, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
     }



    for (float i = -1; i <= 1; i += 0.5) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(i, road_height, (road_length / 2) - 1));
        scale = glm::scale(model, glm::vec3(0.5, 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
     }


    for (float i = -1; i <= 1; i += 0.5) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(i, road_height, -(road_length / 2) + 4));
        scale = glm::scale(model, glm::vec3(0.5, 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
    }





    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, road_height, 5 / 2));
    scale = glm::scale(model, glm::vec3(1.5, 0.3, road_length - 1 - 5));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShader);
    //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
    grass->drawCubeWithTexture(lightingShaderWithTexture, model);




    float field_width = 40;
    ////field

    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(road_width / 2 + field_width / 2 + 0.5, 0, 0));
    //scale = glm::scale(model, glm::vec3(field_width, road_height + 0.5, road_length - 2));
    //translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    //model = alTogether * translate2 * scale * translate;
    //useShaderProgram(lightingShader);
    ////drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
    /////field->drawCubeWithTexture(lightingShaderWithTexture, model);
    //grass->drawCubeWithTexture(lightingShaderWithTexture, model);


    //fieldInside

    float field_width1 = 200;
    float height = 5;
    float field_length = 200;

    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(30, 0, -field_length / 2 + road_length / 2));
    //scale = glm::scale(model, glm::vec3(field_width1, -height, field_length));
    //translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    //model = alTogether * translate2 * scale * translate;
    //useShaderProgram(lightingShader);
    //drawCube(cubeVAO, lightingShader, model, 0.6, 0.4, 0.2);
    ////field2->drawCubeWithTexture(lightingShaderWithTexture, model);



    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(-(road_width / 2 + 2 + 0.5), 0, 0));
    //scale = glm::scale(model, glm::vec3(4, road_height + 0.5, road_length - 2));
    //translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    //model = alTogether * translate2 * scale * translate;
    //useShaderProgram(lightingShader);
    ////drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
    ////field1->drawCubeWithTexture(lightingShaderWithTexture, model);
    //grass->drawCubeWithTexture(lightingShaderWithTexture, model);


    for (float i = 0; i <= field_width; i += 0.5) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(i + road_width / 2 + 0.5, 0, (road_length / 2) - 1));
        scale = glm::scale(model, glm::vec3(0.5, road_height + 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
        //grass->drawCubeWithTexture(lightingShaderWithTexture, model);
    }


    for (float i = 0; i <= field_width; i += 0.5) {
        model = glm::mat4(1.0f);
        translate = glm::mat4(1.0f);
        translate2 = glm::mat4(1.0f);
        scale = glm::mat4(1.0f);
        translate2 = glm::translate(model, glm::vec3(-(i + road_width / 2 + 0.5), 0, (road_length / 2) - 1));
        scale = glm::scale(model, glm::vec3(0.5, road_height + 0.5, 1));
        translate = glm::translate(model, glm::vec3(-0.5, 0, 0));
        model = alTogether * translate2 * scale * translate;
        useShaderProgram(lightingShader);
        //drawCube(cubeVAO, lightingShader, model, 0.471, 0.196, 0.039);
        divider->drawCubeWithTexture(lightingShaderWithTexture, model);
    }

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    glm::mat4 rotate = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(1, 1, 1));

}

float angle = 0;
float angle2 = 0;
float var = 0;
void pond(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether) {
    float pondh = 0.02;
    float pondl = 70;
    float pondw = 30;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(pondw, pondh, pondl));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate * scale;
    pondtex->drawCubeWithTexture(lightingShaderWithTexture, model);
    //, 0.0f, 0.5f, 0.7f
    float fishlength = 3;
    float fishwidth = 2;
    float fishheight = 0.01;

    moveFish += sign * 1.5;
    if (moveFish >= 60)
    {
        sign *= -1;
        angle += 180;
    }
    if (moveFish <= 5)
    {
        sign *= -1;
        angle += 180;
    }
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.8, 1.5, 1.5));
    translate = glm::translate(model, glm::vec3(10.0, 0.08, moveFish));
    glm::mat4 rev = glm::rotate(glm::mat4(1.0f), glm::radians(angle), glm::vec3(1.0f, 0.0f, 0.0f));
    glm::mat4 rot = glm::rotate(glm::mat4(1.0f), glm::radians(-270.0f), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rot * rev * scale;
    //drawFish(lightingShader, model);
    fis->draw(lightingShader, model);

    angle2 += 10;
    if (var > 1)
        var -= 0.01;
    else
        var += 0.01;
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.8, 1.5, 1.5));
    translate = glm::translate(model, glm::vec3(20.0, 0.00+var, 30));
    translate2 = glm::translate(model, glm::vec3(0.0, 2, 0));
    rot = glm::rotate(glm::mat4(1.0f), glm::radians(angle2), glm::vec3(1.0f, 0.0f, 0.0f));
    model = alTogether * translate * rot * translate2 * scale;
    //drawFish(lightingShader, model);
    fis->draw(lightingShader, model);

    // Boarder
    float bl = 30;
    float bl2 = 70;
    float bw = 1;
    float bh = 1;
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bl, bw, bh));
    translate = glm::translate(model, glm::vec3(0, 0, 0-bh));
    model = alTogether * translate * scale;
    cube.drawCube2(lightingShader, model, 1.0f, 1.0f, 1.0f);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bh, bw, bl2+bh));
    translate = glm::translate(model, glm::vec3(bl, 0, 0 - bh));
    model = alTogether * translate * scale;
    cube.drawCube2(lightingShader, model, 1.0f, 1.0f, 1.0f);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bl+bw, bw, bh));
    translate = glm::translate(model, glm::vec3(0, 0, bl2));
    model = alTogether * translate * scale;
    cube.drawCube2(lightingShader, model, 1.0f, 1.0f, 1.0f);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bh, bw, bl2 + bh));
    translate = glm::translate(model, glm::vec3(0, 0, 0 - bh));
    model = alTogether * translate * scale;
    cube.drawCube2(lightingShader, model, 1.0f, 1.0f, 1.0f);
}
void chair(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{
    float backheight = 0.5;
    float backwidth = 0.03;
    float backlength = 0.5;

    float hight = 0.04;
    float length = 0.5;
    float width = 0.3;
    float supporthight = 0.5;
    float supportwidth = 0.04;
    float supportlength = 0.04;
    useShaderProgram(lightingShader);

    // leg1
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(supportlength, supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);

    //leg2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(supportlength, supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(length-supportwidth, 0.0, 0.0));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);
    
    //leg3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(supportlength, supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, width-supportlength));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);

    //leg4
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(supportlength, supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(length-supportwidth, 0.0, width - supportlength));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);

    //base
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length, hight, width));
    translate = glm::translate(model, glm::vec3(0.0, supporthight, 0.0));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);

    //back
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(backlength, backheight, backwidth));
    translate = glm::translate(model, glm::vec3(0.0, supporthight, width - supportlength));
    model = alTogether * translate * scale;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);

    
}
void Fan(Cube& cube, Shader& lightingShader, glm::mat4 alTogether)
{
    float bladel = 1.5;
    float bladew = 0.2;
    float bladeh = 0.01;
    
    //glm::mat4 modelForSphere = glm::mat4(1.0f);
    //modelForSphere = glm::translate(alTogether, glm::vec3(1.5f, 1.2f, 0.5f));
    //sphere->setColor(glm::vec3(0.5, 0.2, 0.5));
    //sphere->drawSphere(lightingShader, modelForSphere);

    // Center
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(0.27, 0.3, 0.27));
    translate = glm::translate(model, glm::vec3(-0.67, 0.0, -0.4));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.2, 0.1, 0.5);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.0, 0.0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(45.0f+rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(165.0f+ rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);


    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(bladel, bladeh, bladew));
    translate = glm::translate(model, glm::vec3(0.01, 0.01, 0.0));
    rotateM = glm::rotate(model, glm::radians(285.0f+ rotateFan), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 1.0f, 0.84f, 0.0f);

    ////leg1
    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(0, 0.02, 0));
    //scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    //translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    //model = alTogether * translate2 * scale * translate;
    //cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);
}
void bed(Cube &cube, Shader & lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{
    float baseHeight = 0.3;
    float width = 1;
    float length = 2;
    float pillowWidth = 0.3;
    float pillowLength = 0.15;
    float blanketWidth = 0.8;
    float blanketLength = 0.7;
    float headHeight = 0.6;

    //base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * scale * translate;
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //float r = 0.345f;
    //float g = 0.171f;
    //float b = 0.475f;
    //cube.setMaterialisticProperty(glm::vec3(r, g, b), glm::vec3(r, g, b), glm::vec3(0.5f, 0.5f, 0.5f));
    //cube.drawCubeWithMaterialisticProperty(lightingShader, model);
    
    //foam
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight, 0));
    scale = glm::scale(model, glm::vec3(width, 0.06, length));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);

    //pillow 1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((width / 2) - (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, .3, .1, .8);

    //pillow 2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3((-width / 2) + (0.1 + pillowWidth / 2), baseHeight + 1 * 0.06, (length / 2) - (0.025 + pillowWidth / 2)));
    scale = glm::scale(model, glm::vec3(pillowWidth, 0.04, pillowLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, .3, .1, .8);

    //blanket
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, baseHeight + 1 * 0.06, -(length / 2 - 0.025) + blanketLength / 2));
    scale = glm::scale(model, glm::vec3(blanketWidth, 0.015, blanketLength));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.541, 0.169, 0.886);

    //head
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, (length / 2 - 0.02 / 2) + 0.02));
    scale = glm::scale(model, glm::vec3(width, headHeight, 0.02));
    translate = glm::translate(model, glm::vec3(-0.5, 0, -0.5));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.545, 0.271, 0.075);
    

}
void table(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{

    float hight = 0.04;
    float length = 0.5;
    float width = 0.3;
    float supporthight = 0.7;
    float supportwidth = 0.04;
    float supportlength = 0.04;

    // base
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(length, hight, length));
    translate = glm::translate(model, glm::vec3(0.0, supporthight, 0.0));
    model = alTogether * scale * translate;
    //drawCube(cubeVAO, lightingShader, model, 0.545, 0.271, 0.075);
    useShaderProgram(lightingShader);
    cube.drawCube2(lightingShader, model, 0.345, 0.171, 0.475);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);


    //leg1
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0.02, 0));
    scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);

    //leg2
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0.02, length - 0.05));
    scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);

    //leg3
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(length - supportwidth, 0.02, length - 0.05));
    scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);

    //leg4
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(length - supportwidth, 0.02, supportlength));
    scale = glm::scale(model, glm::vec3(supportlength, -supporthight, supportwidth));
    translate = glm::translate(model, glm::vec3(0.0, 0.0, 0.0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCube2(lightingShader, model, 0.804, 0.361, 0.361);

}
void door(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{
    float dy = 3.0;
    float dx = 1.4;
    float dz = 0.2;
    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    glm::mat4 rotate = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(dx, dy, dz));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    glm::mat4 rotateM = glm::rotate(model, glm::radians(doorangle), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * rotateM * scale * translate;
    roomdoor->drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.4,0.5,0.1);
}
void stair(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{
    float step = 0.02;
    float h = 0.02;
    //cube.drawCube2(lightingShader, alTogether, 1, 1, 1);
}
void room(Cube& cube, Shader& lightingShaderWithTexture, Shader& lightingShader, glm::mat4 alTogether)
{
    float wx = 2;
    float wy = 2;
    float wz = 0.15;

    /// floor
    float baseHeight = 0.1;
    float width = 6.0;
    float length = 6.5;
    float balconyz = 3;

    float balconyheight = wy / 2 - 0.05;

    glm::mat4 model = glm::mat4(1.0f);
    glm::mat4 translate = glm::mat4(1.0f);
    glm::mat4 translate2 = glm::mat4(1.0f);
    glm::mat4 scale = glm::mat4(1.0f);
    scale = glm::scale(model, glm::vec3(width, baseHeight, length + balconyz));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    roomfloor->drawCubeWithTexture(lightingShaderWithTexture, model);    
    useShaderProgram(lightingShader);
    //cube.drawCube2(lightingShader, model, 0.3, 0.7, 1);

    // Staircase
    //model = glm::mat4(1.0f);
    //translate = glm::mat4(1.0f);
    //translate2 = glm::mat4(1.0f);
    //scale = glm::mat4(1.0f);
    //translate2 = glm::translate(model, glm::vec3(-0.02, 0, length));
    //scale = glm::scale(model, glm::vec3(baseHeight, width, baseHeight));
    //translate = glm::translate(model, glm::vec3(0, 0, 0));
    //model = alTogether * translate2 * scale * translate;
    //useShaderProgram(lightingShaderWithTexture);
    //cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    ////cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);

    // balcony height
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, length + balconyz));
    scale = glm::scale(model, glm::vec3(baseHeight, width, baseHeight));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    
    // balcony support
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, length + balconyz));
    scale = glm::scale(model, glm::vec3(length, balconyheight, baseHeight));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    walltex->drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    //useShaderProgram(lightingShader);

    
    // door
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, length-0.2));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    door(cube, lightingShaderWithTexture, lightingShader, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    //useShaderProgram(lightingShader);

    // doorwall
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 3.0, length));
    scale = glm::scale(model, glm::vec3(width, width-3.0, baseHeight));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    useShaderProgram(lightingShader);

    // doorwall
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(1.4, 0.0, length));
    scale = glm::scale(model, glm::vec3(width-1.4, length - 3.0, baseHeight));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    useShaderProgram(lightingShader);
    
    // Roof
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, width-0.1, balconyz));
    scale = glm::scale(model, glm::vec3(width, baseHeight+0.1, length));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.5, 0.361, 0.361);
    useShaderProgram(lightingShader);

    // Window
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(width/2-wx/2, wy/2, 0-0.05));
    scale = glm::scale(model, glm::vec3(wx, wy, wz+0.1));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    roomwindow->drawCubeWithTexture(lightingShaderWithTexture, model);
    //useShaderProgram(lightingShader);
    //cube.drawCube2(lightingShader, model, 0.5, 0.0, 0.361);

    // Window
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(width / 2 - wx / 2, wy / 2, length - 0.05));
    scale = glm::scale(model, glm::vec3(wx, wy, wz + 0.1));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    roomwindow->drawCubeWithTexture(lightingShaderWithTexture, model);
    //useShaderProgram(lightingShader);
    //cube.drawCube2(lightingShader, model, 0.5, 0.0, 0.361);

    /// right wall
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(width, 0, 0));
    scale = glm::scale(model, glm::vec3(baseHeight, width, length));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);


    /// left wall
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, 0));
    scale = glm::scale(model, glm::vec3(baseHeight, width, length));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);


    /// back wall
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate2 = glm::mat4(1.0f);
    scale = glm::mat4(1.0f);
    translate2 = glm::translate(model, glm::vec3(0, 0, 0));
    scale = glm::scale(model, glm::vec3(length, width, baseHeight));
    translate = glm::translate(model, glm::vec3(0, 0, 0));
    model = alTogether * translate2 * scale * translate;
    useShaderProgram(lightingShaderWithTexture);
    cube.drawCubeWithTexture(lightingShaderWithTexture, model);
    //cube.drawCube2(lightingShader, model, 0.3, 0.7, 1);
    useShaderProgram(lightingShader);
    

    float up = 0.03;
    // Bed
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(1.2, up, 2.0));
    glm::mat4 rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    bed(cube, lightingShaderWithTexture, lightingShader, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length - 1.6, up, 2.0));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    bed(cube, lightingShaderWithTexture, lightingShader, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(1.2, up, width - 1.5));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    bed(cube, lightingShaderWithTexture, lightingShader, model);

    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length - 1.6, up, width - 1.5));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    bed(cube, lightingShaderWithTexture, lightingShader, model);
    float sam = 0.08;
    // table
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length / 2 - width / 6, 4*up+0.5 + sam, 0.1));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    table(cube, lightingShaderWithTexture, lightingShader, model);

    //chair
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length / 2 - width / 6, up+0.0, 0.1 + 1));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    chair(cube, lightingShaderWithTexture, lightingShader, model);


    // table
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length / 2 + width / 8, 4*up + 0.5+sam, 0.1));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    table(cube, lightingShaderWithTexture, lightingShader, model);

    //chair
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length / 2 + width / 8, up + 0.0, 0.1 + 1));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    chair(cube, lightingShaderWithTexture, lightingShader, model);
    // table
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(0.1, 0.5+sam+ 4*up, width / 2));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    table(cube, lightingShaderWithTexture, lightingShader, model);

    //chair
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(0.1+0.7, 0.0, up + width / 2 + 0.5));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    chair(cube, lightingShaderWithTexture, lightingShader, model);

    //table
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length - 1.0, 0.5+sam + 4*up, width / 2));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    table(cube, lightingShaderWithTexture, lightingShader, model);

    //chair
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(length - 0.1 - 1.0, up + 0.0, width / 2));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(-90.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    chair(cube, lightingShaderWithTexture, lightingShader, model);

    // Fan
    model = glm::mat4(1.0f);
    translate = glm::mat4(1.0f);
    translate = glm::translate(model, glm::vec3(1.2+2, width-0.2, 2.0));
    rotateYMatrix = glm::rotate(glm::mat4(1.0f), glm::radians(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    model = alTogether * translate * rotateYMatrix * glm::mat4(1.0f);
    Fan(cube, lightingShader, model);
}


// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
        camera.ProcessKeyboard(FORWARD, deltaTime + 0.4);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
        camera.ProcessKeyboard(BACKWARD, deltaTime + 0.4);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
        camera.ProcessKeyboard(LEFT, deltaTime + 0.4);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
        camera.ProcessKeyboard(RIGHT, deltaTime + 0.4);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS) {
        camera.ProcessKeyboard(UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS) {
        camera.ProcessKeyboard(DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_UP, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS) {
        camera.ProcessKeyboard(P_DOWN, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS) {
        camera.ProcessKeyboard(Y_RIGHT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_LEFT, deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) {
        camera.ProcessKeyboard(R_RIGHT, deltaTime);
    }

    /*if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
    {
        if (rotateAxis_X) rotateAngle_X -= 0.1;
        else if (rotateAxis_Y) rotateAngle_Y -= 0.1;
        else rotateAngle_Z -= 0.1;
    }*/
    //if (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) translate_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) translate_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS) translate_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_J) == GLFW_PRESS) translate_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) translate_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) translate_Z -= 0.001;
    ////if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) scale_X += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS) scale_X -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS) scale_Y += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS) scale_Y -= 0.001;
    //if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) scale_Z += 0.001;
    //if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) scale_Z -= 0.001;

    if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS)
    {
        rotateAngle_X += 0.1;
        rotateAxis_X = 1.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS)
    {
        rotateAngle_Y += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 1.0;
        rotateAxis_Z = 0.0;
    }
    if (glfwGetKey(window, GLFW_KEY_Z) == GLFW_PRESS)
    {
        rotateAngle_Z += 0.1;
        rotateAxis_X = 0.0;
        rotateAxis_Y = 0.0;
        rotateAxis_Z = 1.0;
    }

    /*if (glfwGetKey(window, GLFW_KEY_H) == GLFW_PRESS)
    {
        eyeX += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
    {
        eyeX -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS)
    {
        eyeZ += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        eyeZ -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS)
    {
        eyeY += 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }
    if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
    {
        eyeY -= 2.5 * deltaTime;
        basic_camera.changeEye(eyeX, eyeY, eyeZ);
    }*/
    if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS)
    {
        pointlight2.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_V) == GLFW_PRESS)
    {
        pointlight2.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_B) == GLFW_PRESS)
    {
        pointlight3.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_N) == GLFW_PRESS)
    {
        pointlight3.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS)
    {
        pointlight1.turnOn();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS)
    {
        pointlight1.turnOff();
        // pointlight3.turnOff();
        // pointlight4.turnOff();

    }
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnAmbientOn();
        if (pointlight2.isOn())
            pointlight2.turnAmbientOn();
        if (pointlight3.isOn())
            pointlight3.turnAmbientOn();
        //pointlight4.turnDiffuseOn();
        //diffuseToggle = !diffuseToggle;
    //}
    }
    if (glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnAmbientOff();
        if (pointlight2.isOn())
            pointlight2.turnAmbientOff();
        if (pointlight3.isOn())
            pointlight3.turnAmbientOff();
        //pointlight4.turnDiffuseOff();
        //diffuseToggle = !diffuseToggle;
    //}
    }

    if (glfwGetKey(window, GLFW_KEY_3) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnDiffuseOn();
        if (pointlight2.isOn())
            pointlight2.turnDiffuseOn();
        if (pointlight3.isOn())
            pointlight3.turnDiffuseOn();
        //pointlight4.turnAmbientOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_4) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnDiffuseOff();
        if (pointlight2.isOn())
            pointlight2.turnDiffuseOff();
        if (pointlight3.isOn())
            pointlight3.turnDiffuseOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }


    if (glfwGetKey(window, GLFW_KEY_5) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOn();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOn();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOn();
        //pointlight4.turnSpecularOn();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_6) == GLFW_PRESS)
    {
        /*if (diffuseToggle)
        {*/
        /*cout << "1 " << pointlight1.isOn() << endl;
        cout << pointlight2.isOn() << endl;
        cout << pointlight3.isOn() << endl;*/
        if (pointlight1.isOn())
            pointlight1.turnSpecularOff();
        if (pointlight2.isOn())
            pointlight2.turnSpecularOff();
        if (pointlight3.isOn())
            pointlight3.turnSpecularOff();
        //pointlight4.turnSpecularOff();
        //diffuseToggle = !diffuseToggle;
        //}
    }
    if (glfwGetKey(window, GLFW_KEY_G) == GLFW_PRESS)
    {
        rotateFan += 10.0f;
    }
    if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS)
    {
        doorangle += 10;
        doorangle = min(doorangle, 80.0f);
        //std::cout << doorangle << endl;
    }
    if (glfwGetKey(window, GLFW_KEY_L) == GLFW_PRESS)
    {
        doorangle -= 10;
        doorangle = max(doorangle, 0.0f);
    }
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xposIn, double yposIn)
{
    float xpos = static_cast<float>(xposIn);
    float ypos = static_cast<float>(yposIn);

    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(static_cast<float>(yoffset));
}

unsigned int loadTexture(char const* path, GLenum textureWrappingModeS, GLenum textureWrappingModeT, GLenum textureFilteringModeMin, GLenum textureFilteringModeMax)
{
    unsigned int textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char* data = stbi_load(path, &width, &height, &nrComponents, 0);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, textureWrappingModeS);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, textureWrappingModeT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, textureFilteringModeMin);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, textureFilteringModeMax);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    return textureID;
}
