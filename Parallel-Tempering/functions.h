#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

void write_positions(ofstream& file, int amount, double * positions){
// Writes positions to a file (xyz compatible)
  file << amount << std::endl;
  file << "Particle positions" << std::endl;
  for (int i=0; i<amount; ++i){
    file << "Ar" << "\t" << positions[3*i] << "\t" << positions[3*i+1] << "\t" << positions[3*i+2] << std::endl;
  }
}

void print_init(int amount, int refresh, int timesteps, int threads, int threads_copy, int dump, int copies, int reset_steps, int exchange_rate, double dt, double temp, double mass, double friction, double density, double cutoff_radius, double seed, double difference, double second_well, double height, double temp_spacing, bool periodic, bool reset_origin, bool parallel, bool init_pos, bool init_vel, bool double_well, bool out_pos, const char* init_pos_file, const char* init_vel_file, const char* out_pos_file){
// Prints out initial data
  double region_length = pow((amount*mass)/density,1./3.);	
  std::cerr << "# particles:\t\t" << amount << std::endl;
  std::cerr << "refresh rate:\t\t" << refresh << std::endl;
  std::cerr << "temperature:\t\t"<< temp << std::endl;
  std::cerr << "density:\t\t" << density << std::endl;
  std::cerr << "box sidelength:\t\t"<< region_length << std::endl;
  std::cerr << "cutoff radius:\t\t"<< cutoff_radius << std::endl;
  std::cerr << "mass:\t\t\t" << mass << std::endl;
  std::cerr << "friction:\t\t" << friction << std::endl;
  std::cerr << "# timesteps:\t\t" << timesteps << std::endl;
  std::cerr << "timestep:\t\t" << dt << std::endl;
  std::cerr << "# threads:\t\t" << threads << std::endl;
  std::cerr << std::endl;
  std::cerr << "periodic boundary conditions:\t" << std::boolalpha << periodic << std::endl;
  std::cerr << "reset center of mass to origin:\t" << std::boolalpha << reset_origin << std::endl;
  if (reset_origin){
    std::cerr << "  reset rate:\t\t" << reset_steps << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "parallel tempering:\t\t" << std::boolalpha << parallel << std::endl;
  if (parallel){
    std::cerr << "  amount of copies:\t\t" << copies << std::endl;
    std::cerr << "  exchange rate:\t\t" << exchange_rate << std::endl;
    std::cerr << "  threads per copy:\t\t" << threads_copy << std::endl;
    std::cerr << "  temperature spacing factor:\t" << temp_spacing << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "initial position file used:\t" << std::boolalpha << init_pos << std::endl;
  if (init_pos){
    std::cerr << "  initial position file name:\t" << init_pos_file << std::endl;
  }
  std::cerr << "initial velocity file used:\t" << std::boolalpha << init_vel << std::endl;
  if (init_pos){
    std::cerr << "  initial velocity file name:\t" << init_vel_file << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "external double well potential:\t" << std::boolalpha << double_well << std::endl;
  if (double_well){
    std::cerr << "  energy difference between wells:\t" << difference << std::endl;
    std::cerr << "  location second well:\t\t\t" << second_well << std::endl;
    std::cerr << "  height between wells:\t\t\t" << height << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "# steps after values are stored:\t" << dump << std::endl;
}
