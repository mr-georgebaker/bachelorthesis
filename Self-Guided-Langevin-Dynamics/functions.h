#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#define PI 3.14159265

void write_positions(ofstream& file, int amount, double * positions){
// Writes positions to a file (xyz compatible)
  file << amount << std::endl;
  file << "Particle positions" << std::endl;
  for (int i=0; i<amount; ++i){
    file << "Ar" << "\t" << positions[3*i] << "\t" << positions[3*i+1] << "\t" << positions[3*i+2] << std::endl;
  }
}

void print_init(int amount, int refresh, int timesteps, int threads, int threads_copy, int dump, int copies, int reset_steps, int exchange_rate, double dt, double temp, double mass, double friction, double density_in, double region_length, double cutoff_radius, int seed, double difference, double second_well, double height, double temp_spacing, double bond_length, double max_length, double stiffness, double guiding_factor, double average_time, bool periodic, bool periodic_xz, bool reset_origin, bool parallel, bool init_pos, bool init_vel, bool double_well, bool out_pos, bool polymer, bool harmonic, bool fene, bool self_guided, const char* init_pos_file, const char* init_vel_file, const char* out_pos_file){
// Prints out initial data
  double density = density_in;
  if (!polymer && density==0){
    density = 0.6;
  }
  if (region_length==0 || density_in!=0) region_length = pow((amount*mass)/density,1./3.);
  if (amount==0) amount = 512;
  if (refresh==0) refresh = 1;
  if (timesteps==0) timesteps = 1;
  if (threads==0) threads = 1;
  if (threads_copy==0) threads_copy = 1;
  if (dump==0) dump = 1;
  if (reset_steps==0) reset_steps = 1;
  if (dt==0) dt = 0.001;
  if (temp==0) temp = 1;
  if (mass==0) mass = 1;
  if (cutoff_radius==0) cutoff_radius = (pow((amount*mass)/density,1./3.))/2.;
  if (polymer) cutoff_radius = 1.122462;
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
  std::cerr << "self guided langevin dynamics:\t\t" << std::boolalpha << self_guided << std::endl;
  if (self_guided){
    std::cerr << "  guiding factor:\t\t" << guiding_factor << std::endl;
    std::cerr << "  local average time:\t\t" << average_time << std::endl;
  }
  std::cerr << "periodic boundary conditions:\t\t" << std::boolalpha << periodic << std::endl;
  std::cerr << "periodic boundary conditions xz:\t" << std::boolalpha << periodic_xz << std::endl;
  std::cerr << "reset center of mass to origin:\t\t" << std::boolalpha << reset_origin << std::endl;
  if (reset_origin){
    std::cerr << "  reset rate:\t\t" << reset_steps << std::endl;
  }
  std::cerr << std::endl;
  std::cerr << "polymer:\t\t\t" << std::boolalpha << polymer << std::endl;
  if (polymer){
	std::cerr << "  harmonic potential:\t\t" << std::boolalpha << harmonic << std::endl;
	if (harmonic){
      std::cerr << "    mean bond length:\t\t" << bond_length << std::endl;
    }
    std::cerr << "  FENE potential:\t\t" << std::boolalpha << fene << std::endl;
    if (fene){
	  std::cerr << "    max bond length:\t\t" << max_length << std::endl;
	  std::cerr << "    spring stiffness:\t\t" << stiffness << std::endl;
    }
    std::cerr << "  relaxation time:\t\t" << (friction*(amount-1)*(amount-1)*bond_length*bond_length)/(3*PI*PI*temp) << std::endl;
    std::cerr << "  steps needed for relaxation:\t" << (friction*(amount-1)*(amount-1)*bond_length*bond_length)/(3*PI*PI*temp*dt) << std::endl;
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
