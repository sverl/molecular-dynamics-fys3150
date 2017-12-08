#include "lennardjones.h"
#include "system.h"
#include <math.h>

double LennardJones::potentialEnergy(System& system) {
  m_potentialEnergy = 0;

  int n_atoms = system.atoms().size();
  vec3 cell_size_half = .5 * system.systemSize();
  vec3 cell_size = system.systemSize();

  for (int i = 0; i < n_atoms; i++) {
    Atom* atom_i = system.atoms()[i];
    for (int j = i + 1; j < n_atoms; j++) {
      Atom* atom_j = system.atoms()[j];

      vec3 r_i = atom_i->image;
      vec3 r_j = atom_j->image;
      vec3 r = r_j- r_i;

      // Minimum image
      for (int k = 0; k < 3; k++) {
        if (r[k] > cell_size_half[k]) {
          r[k] -= cell_size[k];
        } else if (r[k] <= cell_size_half[k]) {
          r[k] += cell_size[k];
        }
      }

      float R = r.length();

      m_potentialEnergy +=
          4 * m_epsilon * (pow(m_sigma / R, 12) - pow(m_sigma / R, 6));
    }
  }

  return m_potentialEnergy;
}

double LennardJones::sigma() const {
  return m_sigma;
}

void LennardJones::setSigma(double sigma) {
  m_sigma = sigma;
}

double LennardJones::epsilon() const {
  return m_epsilon;
}

void LennardJones::setEpsilon(double epsilon) {
  m_epsilon = epsilon;
}

void LennardJones::calculateForces(System& system) {
  int n_atoms = system.atoms().size();
  vec3 cell_size_half = .5 * system.systemSize();
  vec3 cell_size = system.systemSize();

  for (int i = 0; i < n_atoms; i++) {
    Atom* atom_i = system.atoms()[i];
    for (int j = i + 1; j < n_atoms; j++) {
      Atom* atom_j = system.atoms()[j];

      vec3 r_i = atom_i->image;
      vec3 r_j = atom_j->image;
      vec3 r = r_j- r_i;

      // Minimum image
      for (int k = 0; k < 3; k++) {
        if (r[k] > cell_size_half[k]) {
          r[k] -= cell_size[k];
        } else if (r[k] <= -cell_size_half[k]) {
          r[k] += cell_size[k];
        }
      }

      float R = r.length();

      vec3 F = 24 * m_epsilon *
               (2 * pow(m_sigma / R, 12) - pow(m_sigma / R, 6)) * r / (R * R) ;

      atom_i->force -= F;
      atom_j->force += F;
    }
  }

}
