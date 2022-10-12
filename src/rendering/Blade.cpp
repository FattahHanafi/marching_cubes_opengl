#include <glm/glm.hpp>

class Blade {
  public:
    Blade(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
    {
        this->old_vertices[0] = v0;
        this->old_vertices[1] = v1;
        this->old_vertices[2] = v2;
        this->old_vertices[3] = v3;

        this->new_vertices[0] = v0;
        this->new_vertices[1] = v1;
        this->new_vertices[2] = v2;
        this->new_vertices[3] = v3;
    }

    void UpdateLocation(glm::vec3 v0, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3)
    {
        old_vertices[0] = new_vertices[0];
        old_vertices[1] = new_vertices[1];
        old_vertices[2] = new_vertices[2];
        old_vertices[3] = new_vertices[3];

        new_vertices[0] = v0;
        new_vertices[1] = v1;
        new_vertices[2] = v2;
        new_vertices[3] = v3;
    }

    void Transform(glm::vec3 trans)
    {
		old_vertices[0] = new_vertices[0];
		old_vertices[1] = new_vertices[1];
		old_vertices[2] = new_vertices[2];
		old_vertices[3] = new_vertices[3];
        new_vertices[0] = old_vertices[0] + trans;
        new_vertices[1] = old_vertices[1] + trans;
        new_vertices[2] = old_vertices[2] + trans;
        new_vertices[3] = old_vertices[3] + trans;
    }

    glm::vec3 old_vertices[4];
    glm::vec3 new_vertices[4];
};
