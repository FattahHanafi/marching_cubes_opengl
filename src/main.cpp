/**
 * Copyright (C) 2018 Tomasz Ga³aj
 **/

#include <cstdint>
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <chrono>
#include <glm/glm.hpp>
#include <iostream>
#include <vector>

#include "rendering/MarchingCube.hpp"
#include "rendering/Mesh.h"
#include "rendering/Model.h"
#include "rendering/Shader.h"
#include "rendering/Texture.h"

GLFWwindow* window;
const int WINDOW_WIDTH = 1280;
const int WINDOW_HEIGHT = 800;

Model* bucket = nullptr;
Shader* shader = nullptr;
Texture* bucket_texture = nullptr;
Texture* ground_texture = nullptr;

unsigned int VBO, EBO, VAO;

MarchingCube mc(200, 100, 100, 0.01f, 0.01f, 0.01f, -0.5f, -0.25f, -0.5f);
Blade blade(glm::vec3(0.3f, 0.1f, 0.1f), glm::vec3(0.3f, 0.1f, 0.0f), glm::vec3(0.3f, -0.1f, 0.0f), glm::vec3(0.3f, -0.1f, 0.1f));
std::vector<Vertex> ground_vertics;
std::vector<unsigned int> ground_indices;

/* Matrices */
// glm::vec3 cam_position = glm::vec3(0.0f, 1.0f, 1.2f);
// glm::vec3 cam_look_at  = glm::vec3(0.0f, 0.5f, 0.0f);
// glm::vec3 cam_up       = glm::vec3(0.0f, 1.0f, 0.0f);

glm::vec3 cam_position = glm::vec3(-1.5f, -1.5f, +1.5f);
glm::vec3 cam_look_at = glm::vec3(0.0f, 0.0f, 0.0f);
glm::vec3 cam_up = glm::vec3(0.0f, 0.0f, 1.0f);

glm::mat4 world_matrix = glm::mat4(1.0f);
glm::mat4 view_matrix = glm::lookAt(cam_position, cam_look_at, cam_up);
// glm::mat4 projection_matrix = glm::perspectiveFov(glm::radians(30.0f), float(WINDOW_WIDTH), float(WINDOW_HEIGHT), 0.1f, 10.0f);
glm::mat4 projection_matrix = glm::ortho(-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 10.0f);
void window_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
    projection_matrix = glm::perspectiveFov(glm::radians(60.0f), float(width), float(height), 0.1f, 10.0f);

    if (shader != nullptr) {
        shader->setUniformMatrix4fv("viewProj", projection_matrix * view_matrix);
    }
}

int init()
{
    /* Initialize the library */
    if (!glfwInit()) return -1;

    /* Create a windowed mode window and its OpenGL context */
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    window = glfwCreateWindow(WINDOW_WIDTH, WINDOW_HEIGHT, "Marching Cubes", nullptr, nullptr);

    if (!window) {
        glfwTerminate();
        return -1;
    }

    /* Make the window's context current */
    glfwMakeContextCurrent(window);

    glfwSetWindowSizeCallback(window, window_size_callback);

    /* Initialize glad */
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    /* Set the viewport */
    glClearColor(0.6784f, 0.8f, 1.0f, 1.0f);
    glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

    glEnable(GL_DEPTH_TEST);

    return true;
}

int loadContent()
{
    bucket = new Model("res/models/model.obj");

	//for(uint32_t idx = 0; idx < mc.getVerticesCount(); ++idx)
	//	mc.vertices[idx] = 1;

	for (int32_t i = 1; i <= 100; ++i) 
        for (int32_t j = 1; j <= 50; ++j)
			for (int32_t k = 0; k <= 50; ++k)
			{
				mc.setVertex(i, j, k, k < 50);
			}

	mc.UpdateCubes();
	mc.CreateModel(&ground_vertics, &ground_indices);
	
	/* Create and apply basic shader */
    shader = new Shader("Basic.vert", "Basic.frag");
    shader->apply();

    shader->setUniformMatrix4fv("world", world_matrix);
    shader->setUniformMatrix3fv("normalMatrix", glm::inverse(glm::transpose(glm::mat3(world_matrix))));
    shader->setUniformMatrix4fv("viewProj", projection_matrix * view_matrix);

    shader->setUniform3fv("cam_pos", cam_position);

    bucket_texture = new Texture();
    bucket_texture->load("res/models/Bucket.png");
    bucket_texture->bind();

    ground_texture = new Texture();
    ground_texture->load("res/models/Ground.png");
    ground_texture->bind();

    return true;
}

void render(float time, float dt)
{
	glm::vec3 relative_pos = glm::vec3(0.3 * (std::cos(0.025f * time) - std::cos(0.025 * (time - dt))), 0.0f, 0.05 * (std::cos(2 * 0.025f * time) - std::cos(2 * 0.025 * (time - dt))));
	glm::vec3 total_pos = glm::vec3(0.05 * std::cos(2 * 0.025f * time),0.0f, -0.3 * std::cos(0.025f * time));
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    /* Draw our triangle */
    world_matrix = glm::rotate(glm::mat4(1.0f), glm::radians(60.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    shader->setUniformMatrix3fv("normalMatrix", glm::inverse(glm::transpose(glm::mat3(world_matrix))));

    blade.Transform(relative_pos);
    mc.CutSoil(&blade);

    /* Update Ground */
    mc.UpdateCubes();

    mc.CreateModel(&ground_vertics, &ground_indices);

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);
    
	glBindVertexArray(VAO);
    // load data into vertex buffers
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, ground_vertics.size() * sizeof(Vertex), &ground_vertics[0], GL_DYNAMIC_DRAW);
    // A great thing about structs is that their memory layout is sequential for all its items.
    // The effect is that we can simply pass a pointer to the struct and it translates perfectly to a glm::vec3/2 array which
    // again translates to 3/2 floats which translates to a byte array.
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, ground_indices.size() * sizeof(unsigned int), &ground_indices[0], GL_DYNAMIC_DRAW);

    // set the vertex attribute pointers
    // vertex Positions
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
    // vertex normals
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, Normal));
    // vertex texture coords
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, TexCoords));
    
	bucket_texture->bind();
    shader->setUniformMatrix4fv("world", glm::translate(glm::rotate(glm::mat4(1.0f), glm::radians(90.0f),glm::vec3(0.0f,1.0f,0.0f)), -total_pos + glm::vec3(-0.05f, 0.0f, 0.0f))); 
    shader->apply();

    bucket->Draw();

    ground_texture->bind();
    shader->setUniformMatrix4fv("world", glm::mat4(1.0f));
    shader->apply();

    glBindVertexArray(VAO);
    glDrawElements(GL_TRIANGLES, ground_indices.size(), GL_UNSIGNED_INT, 0);
   // glBindVertexArray(0);
}

void update()
{
    float startTime = static_cast<float>(glfwGetTime());
	float totalTime = 0;
    /* Loop until the user closes the window */
    while (!glfwWindowShouldClose(window)) {
        /* Update game time value */
        float newTime = static_cast<float>(glfwGetTime());
        float deltaTime = newTime - startTime;
		totalTime += deltaTime;
        /* Render here */
        render(0.01f * totalTime, 0.01f * deltaTime);

        /* Swap front and back buffers */
        glfwSwapBuffers(window);

        /* Poll for and process events */
        glfwPollEvents();
    }
}

int main(void)
{
    if (!init()) return -1;

    if (!loadContent()) return -1;

    update();

    glfwTerminate();

    delete bucket;
    delete shader;
    delete bucket_texture;
    delete ground_texture;

    return 0;
}
