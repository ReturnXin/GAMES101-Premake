#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        for (int i = 0; i < num_nodes; i++)
        {
            Vector2D pos = start + (end - start) / (double)num_nodes * i;
            masses.push_back(new Mass(pos, node_mass, false));
        }

        for (int i = 0; i < num_nodes - 1; i++)
        {
            springs.push_back(new Spring(masses[i], masses[i + 1], k));
        }
//        Comment-in this part when you implement the constructor
        for (auto &i : pinned_nodes) {
            masses[i]->pinned = true;
        }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D a = s->m1->position;
            Vector2D b = s->m2->position;

			float len = (b - a).norm();
            s->m1->forces += -s->k * (a - b) / len * (len - s->rest_length);
			s->m2->forces += -s->k * (b - a) / len * (len - s->rest_length);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                float Kd = 0.005;
                auto a = m->forces / m->mass + gravity - Kd * m->velocity / m->mass;
                m->velocity += a * delta_t;
                m->position += m->velocity * delta_t;

                // TODO (Part 2): Add global damping
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            Vector2D a = s->m1->position;
            Vector2D b = s->m2->position;

            float len = (b - a).norm();
            s->m1->forces += -s->k * (a - b) / len * (len - s->rest_length);
            s->m2->forces += -s->k * (b - a) / len * (len - s->rest_length);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                double damping_factor = 0.00005;
                Vector2D temp_position = m->position;
                // TODO (Part 3.1): Set the new position of the rope mass
                Vector2D a = m->forces / m->mass + gravity;
				m->position = m->position + (1 - damping_factor) * (m->position - m->last_position) + a * delta_t * delta_t;
                m->last_position = temp_position;
                // TODO (Part 4): Add global Verlet damping
            }
            m->forces = Vector2D(0, 0);
        }
    }
}
