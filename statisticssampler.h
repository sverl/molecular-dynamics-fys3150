#ifndef STATISTICSSAMPLER_H
#define STATISTICSSAMPLER_H
#include <fstream>
#include "math/vec3.h"

class System;  // Promise the compiler that this is a class even though we
               // haven't included system.h here

class StatisticsSampler {
 private:
  std::ofstream m_file;
  double m_kineticEnergy = 0;
  double m_potentialEnergy = 0;
  double m_temperature = 0;
  double m_density = 0;
  vec3 m_momentum = vec3(0,0,0);
  float m_meanSquareDev = 0;

 public:
  StatisticsSampler();
  void saveToFile(System& system, string filename);
  void sample(System& system);
  void sampleKineticEnergy(System& system);
  void samplePotentialEnergy(System& system);
  void sampleTemperature(System& system);
  void sampleDensity(System& system);
  void sampleMomentum(System& system);
  void sampleMeanSquareDev(System& system);
  double E_kin() { return m_kineticEnergy; }
  double E_pot() { return m_potentialEnergy; }
  double E_tot() { return m_kineticEnergy + m_potentialEnergy; }
  double temp() { return m_temperature; }
  double density() { return m_density; }
  vec3 momentum(){return m_momentum;}
  double meanSquareDev() { return m_meanSquareDev; }
};
#endif
