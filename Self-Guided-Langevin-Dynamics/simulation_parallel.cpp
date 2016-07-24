/*
* Fluid simulation in a symmetric double well potential
* using the Lennard-Jones-potential and the Langevin-equation.
* A second order integrator and parallel tempering is used.
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <libconfig.h++>
#include "system.h"
#include "functions.h"

using namespace std;
using namespace libconfig;

#ifdef _OPENMP
  #include <omp.h>
#endif

int main(int argc, char* argv[]){

  int amount;
  int timesteps;
  int refresh;
  int threads;
  int threads_copy;
  int dump;
  int copies;
  int reset_steps;
  int exchange_rate;
  int seed;
  double dt;
  double temp;
  double mass;
  double friction;
  double density;
  double region_length;
  double cutoff_radius;
  double difference;
  double second_well;
  double height;
  double temp_spacing;
  double bond_length;
  double max_length;
  double stiffness;
  double* temperatures;
  double guiding_factor;
  double average_time;
  bool self_guided;
  bool periodic;
  bool periodic_xz;
  bool reset_origin;
  bool parallel_tempering;
  bool exchange_temperature;
  bool exchange_momenta;
  bool exchange_coordinates;
  bool init_pos;
  bool init_vel;
  bool double_well;
  bool out_pos;
  bool out_energy;
  bool polymer;
  bool harmonic;
  bool fene;
  const char* init_pos_file;
  const char* init_vel_file;
  const char* out_pos_file;
  const char* out_energy_file;
  System* systems;

  ofstream output_positions;

  /*
  * READ DATA FROM INPUT FILE
  * USING LIBCONFIG 1.5
  */

  Config cfg;
  try{
    cfg.readFile(argv[1]);
  }
  catch(const FileIOException &fioex){
    std::cerr << "I/O error while reading file." << std::endl;
    return(EXIT_FAILURE);
  }
  catch(const ParseException &pex){
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    return(EXIT_FAILURE);
  }
 
  try{
    amount = cfg.lookup("amount_of_particles");
    timesteps = cfg.lookup("timesteps");
    refresh = cfg.lookup("neighborlist_refresh_rate");
    threads = cfg.lookup("threads");
    dump = cfg.lookup("steps_before_dumping_values");
    dt = cfg.lookup("timestep");
    temp = cfg.lookup("temperature");
    mass = cfg.lookup("mass");
    friction = cfg.lookup("friction");
    density = cfg.lookup("density");
    region_length = cfg.lookup("simulation_region_sidelength");
    cutoff_radius = cfg.lookup("cutoff_radius");
    seed = cfg.lookup("seed");
    periodic = cfg.lookup("periodic_boundary_conditions");
    periodic_xz = cfg.lookup("periodic_in_xz");
    reset_origin = cfg.lookup("reset_center_of_mass_to_origin");
    reset_steps = cfg.lookup("steps_between_resetting");
    parallel_tempering = cfg.lookup("parallel_tempering");
    copies = cfg.lookup("amount_of_copies");
    exchange_rate = cfg.lookup("exchange_rate");
    if (copies==0) parallel_tempering = false;
    if (!parallel_tempering){
      copies = 1;
      exchange_rate = 1;
    }
    exchange_temperature = cfg.lookup("exchange_temperature");
    exchange_momenta = cfg.lookup("exchange_momenta");
    exchange_coordinates = cfg.lookup("exchange_coordinates");
    threads_copy = cfg.lookup("threads_per_copy");
    temp_spacing = cfg.lookup("temperature_spacing_factor");
    init_pos = cfg.lookup("initial_position_file");
    init_pos_file = cfg.lookup("initial_position_filename");
    init_vel = cfg.lookup("initial_velocity_file");
    init_vel_file = cfg.lookup("initial_velocity_filename");
    double_well = cfg.lookup("external_double_well_potential");
    difference = cfg.lookup("energy_difference_between_wells");
    second_well = cfg.lookup("location_second_well"); 
    height = cfg.lookup("height_between_wells");
    out_pos = cfg.lookup("output_positions");
    out_pos_file = cfg.lookup("output_positions_file");
    out_energy = cfg.lookup("output_potential_energies");
    out_energy_file = cfg.lookup("output_potential_energies_file");
    polymer = cfg.lookup("polymer");
    harmonic = cfg.lookup("harmonic_potential");
    fene = cfg.lookup("fene_potential");
    bond_length = cfg.lookup("mean_bond_length");
    max_length = cfg.lookup("maximum_bond_length");
    stiffness = cfg.lookup("spring_stiffness");
    self_guided = cfg.lookup("self_guided_langevin_dynamics");
    guiding_factor = cfg.lookup("guiding_factor");
    average_time = cfg.lookup("local_average_time");
  }
  catch(const SettingNotFoundException &nfex){
    cerr << "Error while reading configuration file." << endl;
  }

  /*
  * INITIALIZE
  */

  print_init(amount, refresh, timesteps, threads, threads_copy, dump, copies, reset_steps, exchange_rate, dt, temp, mass, friction, density, region_length, cutoff_radius, seed, difference, second_well, height, temp_spacing, bond_length, max_length, stiffness, guiding_factor, average_time, periodic, periodic_xz, reset_origin, parallel_tempering, init_pos, init_vel, double_well, out_pos, polymer, harmonic, fene, self_guided, init_pos_file, init_vel_file, out_pos_file);

  srand(seed);

  temperatures = new double[copies];
  systems = new System[copies];

  temperatures[0] = temp;
  for (int i=1; i<copies; ++i){
    temperatures[i] = temperatures[i-1]*temp_spacing;
  }
  for (int i=0; i<copies; ++i){
    systems[i].init(amount, timesteps, refresh, threads_copy, dump, reset_steps, dt, temperatures[i], mass, friction, density, region_length, cutoff_radius, bond_length, max_length, stiffness, periodic, periodic_xz,reset_origin, init_pos,init_vel,polymer, harmonic, fene, init_pos_file,init_vel_file, double_well, difference, second_well, height, guiding_factor, average_time);
    systems[i].update_neighborlist();
    systems[i].update_forces();
  }

  output_positions.open(out_pos_file);

  /*
  * SIMULATION
  */

  for (int i=0; i<(timesteps/exchange_rate); ++i){
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int j=0; j<copies; ++j){
      for (int k=0; k<exchange_rate; ++k){
        if (j==0 && (k>=dump || i>=dump) && out_pos){
          //write_positions(output_positions,amount,systems[0].particles.positions);
          std::cout << systems[0].particles.positions[1] << std::endl;
        }
        if (!self_guided){
          systems[j].update_rnd();
          systems[j].update_velocities();
          systems[j].update_positions();
          systems[j].update_forces();
          systems[j].update_velocities();
          systems[j].update_neighborlist();
        }
        else{
          systems[j].update_forces();
          systems[j].update_rnd();
          systems[j].update_guiding_force();
          systems[j].update_scaling_parameter();
          systems[j].update_velocities_self_guided();
          systems[j].update_positions_self_guided();
        }
      }
    }
    #pragma omp barrier
    // Parallel tempering
    for (int j=0; j<(copies-1); ++j){
      double temp_j = systems[j].temperature();
      double temp_j1 = systems[j+1].temperature();
      double delta = (systems[j].potential_energy() - systems[j+1].potential_energy())*(1./temp_j - 1./temp_j1);
      if ((double)rand()/RAND_MAX < exp(delta)){
        //std::cout << "Exchange between replica " << j << " and " << j+1 << " at timestep " << (i+1)*exchange_rate << std::endl;
        if (exchange_coordinates){
          double * cor_j = systems[j].particles.positions;
          double * cor_j1 = systems[j+1].particles.positions;
          systems[j].exchange_coordinates(cor_j1);
          systems[j+1].exchange_coordinates(cor_j);
        }
        if (exchange_momenta){
          double * mom_j = systems[j].particles.velocities;
          double * mom_j1 = systems[j+1].particles.velocities;
          systems[j].exchange_momenta(mom_j1);
          systems[j+1].exchange_momenta(mom_j);
        }
        if (exchange_temperature){
          systems[j].exchange_temperature(temp_j1);
          systems[j+1].exchange_temperature(temp_j);
        }
      }
    }
  }

  /*
  * CLOSE FILES
  */

  output_positions.close();

  return 0;
}
