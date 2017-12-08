#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"

System::System() {}

System::~System() {
  for (Atom* atom : m_atoms) {
    delete atom;
  }
  m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions() {
  for (Atom* atom : m_atoms) {
    for (int i = 0; i < 3; i++) {
      float img = atom->image[i];
      float lim = m_systemSize[i];
      atom->image[i] = img - lim * floor(img / lim);
    }
  }
}

void System::removeTotalMomentum() {
  StatisticsSampler sampler;
  sampler.sampleMomentum(*this);
  vec3 avg_momentum = sampler.momentum() / m_atoms.size();

  for (Atom* atom : m_atoms) {
    atom->velocity -= avg_momentum / (atom->mass());
  }
}

void System::createFCCLattice(int n_cells, double b, double temperature) {
  vec3 cell_origin = vec3(0, 0, 0);

  vec3 r_1 = vec3(0, 0, 0);
  vec3 r_2 = vec3(.5, .5, 0);
  vec3 r_3 = vec3(0, .5, .5);
  vec3 r_4 = vec3(.5, 0, .5);

  vector<vec3> basis = {r_1, r_2, r_3, r_4};

  for (int i = 0; i < n_cells; i++) {
    for (int j = 0; j < n_cells; j++) {
      for (int k = 0; k < n_cells; k++) {
        cell_origin = vec3(i, j, k);
        for (vec3 r : basis) {
          Atom* atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));

          atom->position_init = (cell_origin + r) * b;
          atom->position = atom->position_init;
          atom->image = atom->position_init;
          atom->resetVelocityMaxwellian(temperature);
          m_atoms.push_back(atom);
        }
      }
    }
  }

  float cell_length = b * n_cells;
  setSystemSize(vec3(cell_length, cell_length, cell_length));
}

void System::calculateForces() {
  for (Atom* atom : m_atoms) {
    atom->resetForce();
  }
  m_potential.calculateForces(*this);
}

void System::step(double dt) {
  m_integrator.integrate(*this, dt);
  m_steps++;
  m_time += dt;
}
