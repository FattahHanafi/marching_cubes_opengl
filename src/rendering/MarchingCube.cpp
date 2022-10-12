#include "MarchingCube.hpp"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <list>
#include <unordered_map>
#include <vector>

MarchingCube::MarchingCube(const uint32_t x_steps, const uint32_t y_steps, const uint32_t z_steps, const float x_size, const float y_size,
                           const float z_size, const float x_origin, const float y_origin, const float z_origin)
{
    this->x_steps = x_steps;
    this->y_steps = y_steps;
    this->z_steps = z_steps;

    this->x_size = x_size;
    this->y_size = y_size;
    this->z_size = z_size;

    this->x_origin = x_origin;
    this->y_origin = y_origin;
    this->z_origin = z_origin;

    this->total_cubes = x_steps * y_steps * z_steps;
    this->total_vertices = (x_steps + 1) * (y_steps + 1) * (z_steps + 1);

    this->vertices = new bool[total_vertices];
    this->cubes = new uint8_t[total_cubes];
    this->update_flag = new bool[total_cubes];

    for (uint32_t idx = 0; idx < total_vertices; idx++) {
        vertices[idx] = 0;
    }

    for (uint32_t idx = 0; idx < total_cubes; idx++) {
        cubes[idx] = 0;
        update_flag[idx] = 0;
    }
}

void MarchingCube::setVertex(const uint32_t i, const uint32_t j, const uint32_t k, const bool value)
{
    vertices[getVertexIndex(i, j, k)] = value;
    // UpdateCubeList(i, j, k);
}

uint32_t MarchingCube::getVertexIndex(const uint32_t i, const uint32_t j, const uint32_t k) const
{
    uint32_t idx;
    idx = i * (y_steps + 1) * (z_steps + 1);
    idx += j * (z_steps + 1);
    idx += k;
    return idx;
}

uint32_t MarchingCube::getCubeIndex(const uint32_t i, const uint32_t j, const uint32_t k) const
{
    uint32_t idx;
    idx = i * y_steps * z_steps;
    idx += j * z_steps;
    idx += k;
    return idx;
}

void MarchingCube::getCubeIJK(const uint32_t idx, uint32_t* i, uint32_t* j, uint32_t* k)
{
    *i = idx / (y_steps * z_steps);
    *j = (idx / z_steps) % y_steps;
    *k = idx % z_steps;
}

void MarchingCube::UpdateCubes()
{
    for (uint32_t i = 0; i < x_steps; i++)
        for (uint32_t j = 0; j < y_steps; j++)
            for (uint32_t k = 0; k < z_steps; k++) {
                uint32_t idx = getCubeIndex(i, j, k);
                // if(update_flag[idx] == 0)
                //	continue;
                // update_flag[idx] = 0;
                cubes[idx] = 0;
                cubes[idx] |= vertices[getVertexIndex(i + 0, j + 0, k + 0)] << 0;
                cubes[idx] |= vertices[getVertexIndex(i + 1, j + 0, k + 0)] << 1;
                cubes[idx] |= vertices[getVertexIndex(i + 1, j + 1, k + 0)] << 2;
                cubes[idx] |= vertices[getVertexIndex(i + 0, j + 1, k + 0)] << 3;
                cubes[idx] |= vertices[getVertexIndex(i + 0, j + 0, k + 1)] << 4;
                cubes[idx] |= vertices[getVertexIndex(i + 1, j + 0, k + 1)] << 5;
                cubes[idx] |= vertices[getVertexIndex(i + 1, j + 1, k + 1)] << 6;
                cubes[idx] |= vertices[getVertexIndex(i + 0, j + 1, k + 1)] << 7;
            }

    sum = 0;
    for (uint32_t idx = 0; idx < total_vertices; idx++) {
        sum += vertices[idx];
    }
}

void MarchingCube::UpdateCubeList(const uint32_t i, const uint32_t j, const uint32_t k)
{
    std::list<int32_t> X = {-1, 0};
    std::list<int32_t> Y = {-1, 0};
    std::list<int32_t> Z = {-1, 0};

    if (i == 0) X.pop_front();
    if (i == (x_steps + 1)) X.pop_back();

    if (j == 0) Y.pop_front();
    if (j == (y_steps + 1)) Y.pop_back();

    if (k == 0) Z.pop_front();
    if (k == (z_steps + 1)) Z.pop_back();

    for (int32_t u : X)
        for (int32_t v : Y)
            for (int32_t w : Z) update_flag[getCubeIndex(i + u, j + v, k + w)] = true;
}

uint32_t MarchingCube::getCubesCount() const { return total_cubes; }

uint32_t MarchingCube::getVerticesCount() const { return total_vertices; }

void MarchingCube::CreateModel(std::vector<Vertex>* vertices, std::vector<unsigned int>* indices)
{
    vertices->clear();
    indices->clear();
    Vertex vtx;
    glm::vec3 orig(x_origin, y_origin, z_origin);
    vtx.Normal = glm::vec3(0.0f, 0.0f, 1.0f);
    for (uint32_t I = 0; I < x_steps; ++I)
        for (uint32_t J = 0; J < y_steps; ++J)
            for (uint32_t K = 0; K < z_steps; ++K) {
                uint32_t idx = getCubeIndex(I, J, K);
                for (uint32_t i = 0; i < 16; ++i) {
                    switch (triangles[cubes[idx]][i]) {
                        case 0:
                            vtx.Position = orig + glm::vec3((I + 0.5f) * x_size, (J + 0.0f) * y_size, (K + 0.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 1:
                            vtx.Position = orig + glm::vec3((I + 1.0f) * x_size, (J + 0.5f) * y_size, (K + 0.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 2:
                            vtx.Position = orig + glm::vec3((I + 0.5f) * x_size, (J + 1.0f) * y_size, (K + 0.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 3:
                            vtx.Position = orig + glm::vec3((I + 0.0f) * x_size, (J + 0.5f) * y_size, (K + 0.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 4:
                            vtx.Position = orig + glm::vec3((I + 0.5f) * x_size, (J + 0.0f) * y_size, (K + 1.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 5:
                            vtx.Position = orig + glm::vec3((I + 1.0f) * x_size, (J + 0.5f) * y_size, (K + 1.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 6:
                            vtx.Position = orig + glm::vec3((I + 0.5f) * x_size, (J + 1.0f) * y_size, (K + 1.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 7:
                            vtx.Position = orig + glm::vec3((I + 0.0f) * x_size, (J + 0.5f) * y_size, (K + 1.0f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 8:
                            vtx.Position = orig + glm::vec3((I + 0.0f) * x_size, (J + 0.0f) * y_size, (K + 0.5f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 9:
                            vtx.Position = orig + glm::vec3((I + 1.0f) * x_size, (J + 0.0f) * y_size, (K + 0.5f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 10:
                            vtx.Position = orig + glm::vec3((I + 1.0f) * x_size, (J + 1.0f) * y_size, (K + 0.5f) * z_size);
                            vertices->push_back(vtx);
                            break;
                        case 11:
                            vtx.Position = orig + glm::vec3((I + 0.0f) * x_size, (J + 1.0f) * y_size, (K + 0.5f) * z_size);
                            vertices->push_back(vtx);
                            break;
                    }
                }
            }

    for (uint32_t i = 0; i < vertices->size(); i++) {
        indices->push_back(i);
    }

    for (uint32_t i = 0; i < vertices->size(); i += 3) {
        glm::vec3 nor =
            glm::cross((vertices->at(i + 1).Position - vertices->at(i).Position), (vertices->at(i + 2).Position - vertices->at(i).Position));
        vertices->at(i).Normal = nor;
        vertices->at(i + 1).Normal = nor;
        vertices->at(i + 2).Normal = nor;
        vertices->at(i).TexCoords = glm::vec2(0.0f, 0.0f);
        vertices->at(i + 1).TexCoords = glm::vec2(1.0f, 0.0f);
        vertices->at(i + 2).TexCoords = glm::vec2(1.0f, 1.0f);
    }
}

void MarchingCube::CutSoil(Blade* blade)
{
    float x_min = std::min(blade->new_vertices[0].x, blade->new_vertices[1].x);
    x_min = std::min(x_min, blade->new_vertices[2].x);
    x_min = std::min(x_min, blade->new_vertices[3].x);
    x_min = std::min(x_min, blade->old_vertices[0].x);
    x_min = std::min(x_min, blade->old_vertices[1].x);
    x_min = std::min(x_min, blade->old_vertices[2].x);
    x_min = std::min(x_min, blade->old_vertices[3].x);

    float x_max = std::max(blade->new_vertices[0].x, blade->new_vertices[1].x);
    x_max = std::max(x_max, blade->new_vertices[2].x);
    x_max = std::max(x_max, blade->new_vertices[3].x);
    x_max = std::max(x_max, blade->old_vertices[0].x);
    x_max = std::max(x_max, blade->old_vertices[1].x);
    x_max = std::max(x_max, blade->old_vertices[2].x);
    x_max = std::max(x_max, blade->old_vertices[3].x);

    float y_min = std::min(blade->new_vertices[0].y, blade->new_vertices[1].y);
    y_min = std::min(y_min, blade->new_vertices[2].y);
    y_min = std::min(y_min, blade->new_vertices[3].y);
    y_min = std::min(y_min, blade->old_vertices[0].y);
    y_min = std::min(y_min, blade->old_vertices[1].y);
    y_min = std::min(y_min, blade->old_vertices[2].y);
    y_min = std::min(y_min, blade->old_vertices[3].y);

    float y_max = std::max(blade->new_vertices[0].y, blade->new_vertices[1].y);
    y_max = std::max(y_max, blade->new_vertices[2].y);
    y_max = std::max(y_max, blade->new_vertices[3].y);
    y_max = std::max(y_max, blade->old_vertices[0].y);
    y_max = std::max(y_max, blade->old_vertices[1].y);
    y_max = std::max(y_max, blade->old_vertices[2].y);
    y_max = std::max(y_max, blade->old_vertices[3].y);

    float z_min = std::min(blade->new_vertices[0].z, blade->new_vertices[1].z);
    z_min = std::min(z_min, blade->new_vertices[2].z);
    z_min = std::min(z_min, blade->new_vertices[3].z);
    z_min = std::min(z_min, blade->old_vertices[0].z);
    z_min = std::min(z_min, blade->old_vertices[1].z);
    z_min = std::min(z_min, blade->old_vertices[2].z);
    z_min = std::min(z_min, blade->old_vertices[3].z);

    float z_max = std::max(blade->new_vertices[0].z, blade->new_vertices[1].z);
    z_max = std::max(z_max, blade->new_vertices[2].z);
    z_max = std::max(z_max, blade->new_vertices[3].z);
    z_max = std::max(z_max, blade->old_vertices[0].z);
    z_max = std::max(z_max, blade->old_vertices[1].z);
    z_max = std::max(z_max, blade->old_vertices[2].z);
    z_max = std::max(z_max, blade->old_vertices[3].z);

    uint32_t w = 0;
    glm::vec3 P;
    for (uint32_t i = 0; i <= x_steps; i++) {
        P.x = x_origin + i * x_size;
        if ((P.x < x_max) && (P.x > x_min))
            for (uint32_t j = 0; j <= y_steps; j++) {
                P.y = y_origin + j * y_size;
                if ((P.y < y_max) && (P.y > y_min))
                    for (uint32_t k = 0; k <= z_steps; k++) {
                        P.z = z_origin + k * z_size;
                        if ((P.z < z_max) && (P.z > z_min))
                            if (vertices[getVertexIndex(i, j, k)])
                                if (IsInside(blade, &P)) {
                                    w++;
                                    vertices[getVertexIndex(i, j, k)] = 0;
                                }
                    }
            }
    }

    std::cout << "Updated vertices = " << w << std::endl;
}

bool MarchingCube::IsInside(Blade* blade, glm::vec3* P)
{
    bool s0, s1, s2, s3, s4, s5;
    s0 = sideSign(&blade->old_vertices[0], &blade->old_vertices[1], &blade->old_vertices[2], P);
    s1 = sideSign(&blade->new_vertices[0], &blade->new_vertices[1], &blade->new_vertices[2], P);
    if (s0 == s1) return 0;

    s2 = sideSign(&blade->old_vertices[0], &blade->new_vertices[0], &blade->new_vertices[1], P);
    s3 = sideSign(&blade->old_vertices[3], &blade->new_vertices[3], &blade->new_vertices[2], P);
    if (s2 == s3) return 0;

    s4 = sideSign(&blade->old_vertices[0], &blade->new_vertices[0], &blade->new_vertices[3], P);
    s5 = sideSign(&blade->old_vertices[1], &blade->new_vertices[1], &blade->new_vertices[2], P);
    if (s4 == s5) return 0;

    return 1;
}

bool MarchingCube::sideSign(glm::vec3* A, glm::vec3* B, glm::vec3* C, glm::vec3* P)
{
    glm::vec3 B_ = *B - *A;
    glm::vec3 C_ = *C - *A;
    glm::vec3 P_ = *P - *A;

    // Bx*Cy*Pz - Bx*Cz*Py - By*Cx*Pz + By*Cz*Px + Bz*Cx*Py - Bz*Cy*Px
    float det = +B_.x * C_.y * P_.z;
    det += -B_.x * C_.z * P_.y;
    det += -B_.y * C_.x * P_.z;
    det += +B_.y * C_.z * P_.x;
    det += +B_.z * C_.x * P_.y;
    det += -B_.z * C_.y * P_.x;

    return (det > 0);
}

