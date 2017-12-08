#include "system.h"
#include "statisticssampler.h"
#include "lennardjones.h"
#include <iostream>
#include "math/vec3.h"
#include "unitconverter.h"
#include <iomanip>

using std::ofstream;
using std::cout;
using std::endl;
using std::setw;

StatisticsSampler::StatisticsSampler() {}

void StatisticsSampler::saveToFile(System& system) {
  // Save the statistical properties for each timestep for plotting etc.
  // First, open the file if it's not open already
  if (!m_file.good()) {
    m_file.open("statistics.txt", ofstream::out);
    // If it's still not open, something bad happened...
    if (!m_file.good()) {
      cout << "Error, could not open statistics.txt" << endl;
      exit(1);
    }
  }

  m_file << system.steps() << setw(11)
         << UnitConverter::timeToSI(system.time()) * 1e12 << setw(11)
         << UnitConverter::temperatureToSI(temp()) << setw(11)
         << UnitConverter::energyToEv(E_kin()) << setw(11)
         << UnitConverter::energyToEv(E_pot()) << setw(11)
         << UnitConverter::energyToEv(E_tot()) << setw(11)
         << UnitConverter::lengthToAngstroms(meanSquareDev()) << endl;
}

void StatisticsSampler::sample(System& system) {
  // Here you should measure different kinds of statistical properties and save
  // it to a file.
  sampleKineticEnergy(system);
  samplePotentialEnergy(system);
  sampleTemperature(system);
  sampleDensity(system);
  sampleMeanSquareDev(system);
  saveToFile(system);
}

void StatisticsSampler::sampleKineticEnergy(System& system) {
  m_kineticEnergy = 0;
  for (Atom* atom : system.atoms()) {
    m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
  }
}

void StatisticsSampler::samplePotentialEnergy(System& system) {
  m_potentialEnergy = system.potential().potentialEnergy(system);
}

void StatisticsSampler::sampleTemperature(System& system) {
  m_temperature = 2 * m_kineticEnergy / (3 * system.atoms().size());
}

void StatisticsSampler::sampleDensity(System& system) {}

void StatisticsSampler::sampleMomentum(System& system) {
  m_momentum.zeros();
  for (Atom* atom : system.atoms()) {
    m_momentum += atom->velocity * atom->mass();
  }
}

void StatisticsSampler::sampleMeanSquareDev(System& system)
{
    m_meanSquareDev = 0;
    for(Atom *atom : system.atoms()) {
        m_meanSquareDev += (atom->position - atom->position_init).lengthSquared();
    }
}
