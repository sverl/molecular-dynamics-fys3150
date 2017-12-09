#include "math/random.h"
#include "lennardjones.h"
#include "velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>
#include <iomanip>

using namespace std;

int main(int argc, char** argv) {
  int n_cells = 5;
  double temp_init =
      UnitConverter::temperatureFromSI(300.0);          // measured in Kelvin
  double b = UnitConverter::lengthFromAngstroms(5.26);  // measured in angstroms
  int timestep_max = 10000;

  // If a first argument is provided, it is the number of unit cells
  if (argc > 1)
    n_cells = atoi(argv[1]);
  // If a second argument is provided, it is the initial temperature (measured
  // in kelvin)
  if (argc > 2)
    temp_init = UnitConverter::temperatureFromSI(atof(argv[2]));
  // If a third argument is provided, it is the lattice constant determining the
  // density (measured in angstroms)
  if (argc > 3)
    b = UnitConverter::lengthFromAngstroms(atof(argv[3]));
  if (argc > 4)
    timestep_max = atoi(argv[4]);

  double dt = UnitConverter::timeFromSI(1e-15);  // Measured in seconds.

  cout << "One unit of length is " << UnitConverter::lengthToSI(1.0)
       << " meters" << endl;
  cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0)
       << " meters/second" << endl;
  cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds"
       << endl;
  cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg"
       << endl;
  cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0)
       << " K" << endl;

  string filename = "statistics_" + to_string(temp_init) + ".txt";

  System system;
  system.createFCCLattice(n_cells, b, temp_init);
  system.potential().setEpsilon(1);
  system.potential().setSigma(UnitConverter::lengthFromAngstroms(3.405));

  system.removeTotalMomentum();

  cout << "Cell size is " << system.systemSize() << endl;

  StatisticsSampler sampler;
  IO movie("movie.xyz");  // To write the state to file

  cout << setw(16) << "Timestep" << setw(16) << "Time" << setw(16)
       << "Temperature" << setw(16) << "KineticEnergy" << setw(16)
       << "PotentialEnergy" << setw(16) << "TotalEnergy" << endl;
  for (int timestep = 0; timestep < timestep_max; timestep++) {
    system.step(dt);
    sampler.sample(system);
    if (timestep % 100 == 0) {
      // Print the timestep every 100 timesteps
      cout << setw(16) << system.steps()   //
           << setw(16) << system.time()    //
           << setw(16) << sampler.temp()   //
           << setw(16) << sampler.E_kin()  //
           << setw(16) << sampler.E_pot()  //
           << setw(16) << sampler.E_tot() << endl;

      movie.saveState(system);
      sampler.saveToFile(system, filename);
    }
  }

  movie.close();

  return 0;
}
