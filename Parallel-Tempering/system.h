#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#define PI 3.14159265

using namespace std;

double distance_periodic(double * xi, double * xj, double * dr, double region_length){
// Returns the squared distance between two points in three dimensions
// for periodic boundary conditions in x,y and z
  double d=0, c, r=region_length/2.;
  c = xi[0]-xj[0];
  if (c>r){
    c -= region_length;
  }
  else if (c<-r){
    c += region_length;
  }
  d += c*c;
  dr[0] = c;
  c = xi[1]-xj[1];
  if (c>r){
    c -= region_length;
  }
  else if (c<-r){
    c += region_length;
  }
  d += c*c;
  dr[1] = c;
  c = xi[2]-xj[2];
  if (c>r){
    c -= region_length;
  }
  else if (c<-r){
    c += region_length;
  }
  d += c*c;
  dr[2] = c;
  return d;
}
double distance_non_periodic(double * xi, double * xj, double * dr){
// Returns the distance squared between two points in three dimensions
  double d=0, c;
  c = xi[0] - xj[0];
  d += c*c;
  dr[0] = c;
  c = xi[1] - xj[1];
  d += c*c;
  dr[1] = c;
  c = xi[2] - xj[2];
  d += c*c;
  dr[2] = c;
  return d;
}

class Particles{
  public:
    double virial;
    double potential_energy;
    double current_temp;
    double init_temp;
    double sigma;
    double * positions;
    double * velocities;
    double * forces;
    double * rnd;
    int ** neighborlist;
    void init(int amount_in, int threads_in, int refresh_in, int reset_steps_in, double dt_in, double temp_in, double mass_in, double friction_in, double region_length_in, double cutoff_radius_in, double difference_in, double second_well_in, double height_in, bool periodic_in, bool reset_origin_in, bool init_pos_in, bool init_vel_in, bool double_well_in, const char* init_pos_file_in, const char* init_vel_file_in);
    void update_neighborlist(void);
    void update_forces(void);
    void update_positions(void);
    void update_positions_verlet(void);
    void update_velocities(void);
    void update_velocities_verlet(void);
    void calc_rnd(void);
  private:
    int amount;
    int threads;
    int refresh;
    int reset_steps;
    int called;
    double dt;
    double mass;
    double friction;
    double region_length;
    double cutoff_radius;
    double verlet_radius;
    double difference;
    double second_well;
    double height;
    double dvi2;
    double dvi6;
    double s3dt;
    double sdt;
    double si12;
    double si3;
    double potential_energy_rc;
    bool periodic;
    bool reset_origin;
    bool init_pos;
    bool init_vel;
    bool double_well;
    const char* init_pos_file;
    const char* init_vel_file;
    void init_positions(void);
    void init_velocities(void);
    double wrap_positions(double pos_in);
    void reset_positions(void);
    void scale_velocities(void);
};
void Particles::init(int amount_in, int threads_in, int refresh_in, int reset_steps_in, double dt_in, double temp_in, double mass_in, double friction_in, double region_length_in, double cutoff_radius_in, double difference_in, double second_well_in, double height_in, bool periodic_in, bool reset_origin_in, bool init_pos_in, bool init_vel_in, bool double_well_in, const char* init_pos_file_in, const char* init_vel_file_in){
  amount = amount_in;
  threads = threads_in;
  refresh = refresh_in;
  reset_steps = reset_steps_in;
  dt = dt_in;
  sdt = sqrt(dt);
  s3dt = sqrt(dt*dt*dt);
  si12 = 1./sqrt(12.);
  si3 = 1./sqrt(3.);
  init_temp = temp_in;
  current_temp = init_temp;
  mass = mass_in;
  friction = friction_in;
  region_length = region_length_in;
  cutoff_radius = cutoff_radius_in;
  verlet_radius = cutoff_radius + 30*dt;
  dvi2 = 1./(verlet_radius*verlet_radius);
  dvi6 = dvi2*dvi2*dvi2;
  potential_energy_rc = 4.*dvi6*(dvi6-1.);
  difference = difference_in;
  second_well = second_well_in;
  height = height_in;
  periodic = periodic_in;
  reset_origin = reset_origin_in;
  init_pos = init_pos_in;
  init_vel = init_vel_in;
  double_well = double_well_in;
  init_pos_file = init_pos_file_in;
  init_vel_file = init_vel_file_in;
  sigma = sqrt((2*init_temp*friction)/mass);
  positions = new double[amount*3];
  velocities = new double[amount*3];
  forces = new double[amount*3];
  neighborlist = new int*[amount-1];
  for (int i=0; i<amount; ++i){
    neighborlist[i] = new int[amount];
  }
  if (friction!=0){
    rnd = new double[amount*6];
  }
  called = 0;
  init_positions();
  init_velocities();
}
void Particles::update_neighborlist(void){
// Updates the neighborlist
  if (called%refresh==0){
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<(amount-1); ++i){
      int l = 0;
      for (int j=(i+1); j<amount; ++j){
        double pos_i[3], pos_j[3], dr[3];
        double distance_ij;
        for (int k=0; k<3; ++k){
          pos_i[k] = positions[3*i+k];
          pos_j[k] = positions[3*j+k];
        }
        if (periodic){
          distance_ij = distance_periodic(pos_i,pos_j,dr,region_length);
        }
        else{
          distance_ij = distance_non_periodic(pos_i,pos_j,dr);
        }
        if (distance_ij < verlet_radius*verlet_radius){
          neighborlist[i][l+1] = j;
          l++;
        }
      }
      neighborlist[i][0] = l;
    }
    #pragma omp barrier
  }
  called++;
}
void Particles::update_forces(void){
// Updates the force based on the Lennard-Jones-potential
// and the external double well potential
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      forces[3*i+j] = 0;
    }
  }
  #pragma omp barrier
  ///////////////////////////////
  /// LENNARD-JONES-POTENTIAL ///
  ///////////////////////////////
  double pot_energy = 0, vir = 0; 
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for reduction(+:vir,pot_energy)
  #endif
  for (int i=0; i<(amount-1); ++i){
    int neighbors = neighborlist[i][0];
    for (int j=1; j<=neighbors; ++j){
      double pos_i[3], pos_j[3], dr[3];
      double d2;
      int index_j = neighborlist[i][j];
      for (int k=0; k<3; ++k){
        pos_i[k] = positions[3*i+k];
        pos_j[k] = positions[3*index_j+k];
      }
      if (periodic){
        d2 = distance_periodic(pos_i,pos_j,dr,region_length);
      }
      else{
        d2 = distance_non_periodic(pos_i,pos_j,dr);
      }
      double di2 = 1./d2;
      double di6 = di2*di2*di2;
      double force = 48.*di6*(di6-0.5)*di2;
      vir += force*sqrt(d2);
      for (int k=0; k<3; ++k){ 
        #pragma omp atomic   
          forces[3*i+k] += force*dr[k];
        #pragma omp atomic
          forces[3*index_j+k] -= force*dr[k];
      }
      pot_energy += 4.*di6*(di6-1.) - potential_energy_rc;
    }
  }
  #pragma omp barrier
  virial = vir;
  potential_energy = pot_energy;
  //////////////////////////////////////
  /// EXTERNAL DOUBLE WELL POTENTIAL ///
  //////////////////////////////////////
  if (double_well){
    double b = height;
    double w = second_well;
    double s = difference;
    double pot_energy = 0, vir = 0; 
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for reduction(+:vir,pot_energy)
    #endif
    for (int i=0; i<amount; ++i){
      double y = positions[3*i+1];
      double force_y = -((2*y*b)/(w*w*w*w))*(y-w)*((y-w)+y)-(s/w);
      forces[3*i+1] += force_y;
      pot_energy += (b/(w*w*w*w))*y*y*(y-w)*(y-w) + ((s*y)/w);
      vir += y*force_y;
    }
    #pragma omp barrier
    virial += vir;
    potential_energy += pot_energy;
  }
} 
void Particles::update_positions(void){
// Updates the positions according to the Langevin integrator;
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double pos = positions[3*i+j] + dt*velocities[3*i+j] + s3dt*sigma*si12*rnd[6*i+j+3];
      if (periodic){
        positions[3*i+j] = wrap_positions(pos);
      }
      else{
        positions[3*i+j] = pos;
      }
    }
  }
  #pragma omp barrier
  if (reset_origin){
    reset_positions();
  }
}
void Particles::update_positions_verlet(void){
// Updates the postions according to the velocity verlet integrator
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double pos = positions[3*i+j] + dt*velocities[3*i+j];
      if (periodic){
        positions[3*i+j] = wrap_positions(pos);
      }
      else{
        positions[3*i+j] = pos;
      }
    }
  }
  #pragma omp barrier
  if (reset_origin){
    reset_positions();
  }
}
void Particles::reset_positions(void){
// Resets the center of mass to the origin 
  if (called%reset_steps==0){
    double center_of_mass = 0;
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for reduction(+:center_of_mass)
    #endif
    for (int i=0; i<amount; ++i){
      center_of_mass += positions[3*i];
    }
    #pragma omp barrier
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<amount; ++i){
      positions[3*i] -= center_of_mass;
    }
  }
}
double Particles::wrap_positions(double pos_in){
// Wraps the position to keep the particle inside the box
  if (pos_in > region_length){
    pos_in -= region_length;
  }
  else if (pos_in < 0){
    pos_in += region_length;
  }
  return pos_in;
}
void Particles::update_velocities(void){
// Updates the velocity according to the Langevin integrator and
// rescales the velocity if the temperature was exchanged
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      velocities[3*i+j] += 0.5*dt*forces[3*i+j]/mass - 0.5*dt*friction*velocities[3*i+j]  + 0.5*sdt*sigma*rnd[6*i+j] - 0.125*dt*dt*friction*(forces[3*i+j]/mass - friction*velocities[3*i+j]) - 0.25*s3dt*friction*sigma*(0.5*rnd[6*i+j] + si3*rnd[6*i+j+3]);
    }
  }
  #pragma omp barrier
  scale_velocities();
}
void Particles::update_velocities_verlet(void){
// Updates the postions according to the velocity verlet integrator
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      velocities[3*i+j] += 0.5*dt*forces[3*i+j]/mass;
    }
  }
  #pragma omp barrier
  scale_velocities();
}
void Particles::scale_velocities(void){
// Scales the velocities to satisfy mean kinetic energy if parallel tempering is used
  if (init_temp/current_temp!=1){
    double scale = sqrt(init_temp/current_temp);
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<amount; ++i){
      for (int j=0; j<3; ++j){
        velocities[3*i+j] *= scale;
      }
    }
    #pragma omp barrier
  }
}
void Particles::calc_rnd(void){
// Calculates the random, normally distributed values needed for the Langevin integrator
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double u = (double)rand()/RAND_MAX;
      double v = (double)rand()/RAND_MAX;
      rnd[6*i+j] = sqrt(-2*log(u))*cos(2*PI*v);
      rnd[6*i+j+3] = sqrt(-2*log(u))*sin(2*PI*v);
    }
  }
}
void Particles::init_positions(void){
// Either reads a file containing the initial positions or sets them on a crystal like grid
  int ix=0, iy=0, iz=0, spacing=1;
  if (!init_pos){
    while ((spacing*spacing*spacing)<amount) spacing++;
    for (int i=0; i<amount; ++i){
      positions[3*i] = ((ix+0.5)*region_length)/spacing;
      positions[3*i+1] = ((iy+0.5)*region_length)/spacing;
      positions[3*i+2] = ((iz+0.5)*region_length)/spacing;
      ix++;
      if (ix==spacing){
        ix = 0;
        iy++;
        if (iy==spacing){
          iy = 0;
          iz++;
        }
      }
    }
  }
  else{
    ifstream init_pos_f (init_pos_file);
    if (init_pos_f.is_open()){
      // Reads from a .xyz file
      // first row: amount of particles
      // second row: description
      // following rows: particle positions
      int amount_file;
      int i=0;
      string description;
      string atom;
      double x,y,z;
      init_pos_f >> amount_file;
      init_pos_f >> description;
      if (amount == amount_file){
        while (init_pos_f >> atom >> x >> y >> z){
          positions[3*i] = x;
          positions[3*i+1] = y;
          positions[3*i+2] = z;
          ++i;
        }
      }
      else{
        std::cerr << "The amount of particles inside the file does not match the specified amount by the configuration file" << std::endl;
        exit(EXIT_FAILURE);
      }
      init_pos_f.close();
    }
    else{
      perror(init_pos_file);
    } 
  }
}
void Particles::init_velocities(void){
// Either reads a file containing the initial velocities or draw them from a Boltzmann-distribution
// and removes net momentum to avoid a center of mass drift
  double mom_x=0, mom_y=0, mom_z=0,u,v,r;
  if (!init_vel){
    for (int i=0; i<amount; ++i){
      for (int j=0; j<3; ++j){
        u = (double)rand()/RAND_MAX;
        v = (double)rand()/RAND_MAX;
        r = sqrt(init_temp/mass)*sqrt(-2*log(u))*cos(2*PI*v);
        velocities[3*i+j] = r;
        if (j==0){
          mom_x += r;
        }
        else if (j==1){
          mom_y += r;
        }
        else{ 
          mom_z += r;
        }
      }
    }
    mom_x /= amount;
    mom_y /= amount;
    mom_z /= amount;
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<amount; ++i){
      velocities[3*i] -= mom_x;
      velocities[3*i+1] -= mom_y;
      velocities[3*i+2] -= mom_z;
    }
    #pragma omp barrier
  }
  else{
    ifstream init_vel_f (init_vel_file);
    if (init_vel_f.is_open()){
      // Reads from a .xyz file
      // first row: amount of particles
      // second row: description
      // following rows: particle velocities
      int amount_file;
      int i=0;
      string description;
      string atom;
      double x,y,z;
      init_vel_f >> amount_file;
      init_vel_f >> description;
      if (amount == amount_file){
        while (init_vel_f >> atom >> x >> y >> z){
          velocities[3*i] = x;
          velocities[3*i+1] = y;
          velocities[3*i+2] = z;
          ++i;
        }
      }
      else{
        std::cerr << "The amount of particles inside the file does not match the specified amount by the configuration file" << std::endl;
        exit(EXIT_FAILURE);
      }
      init_vel_f.close();
    }
    else{
      perror(init_vel_file);
    } 
  }
}

class System{
  public:
    Particles particles;
    void init(int amount_in, int timesteps_in, int refresh_in, int threads_in, int dump_in, int reset_steps_in, double dt_in, double temp_in, double mass_in, double friction_in, double density_in, double cutoff_radius_in, bool periodic_in, bool reset_origin_in, bool init_pos_in, bool init_vel_in, const char* init_pos_file_in, const char* init_vel_file_in, bool double_well_in, double difference_in, double second_well_in, double height_in);
    void update_neighborlist(void);
    void update_forces(void);
    void update_positions(void);
    void update_velocities(void);
    void update_rnd(void);
    double kinetic_energy(void);
    double potential_energy(void);
    double temperature(void);
    double actual_temperature(void);
    double pressure(void);
    void set_temperature(double temp_in);
  private:
    int amount;
    int timesteps;
    int refresh;
    int threads;
    int dump;
    int reset_steps;
    double dt;
    double mass;
    double init_temp;
    double current_temp;
    double friction;
    double density;
    double region_length;
    double cutoff_radius;
    double difference;
    double second_well;
    double height;
    double volume;
    bool periodic;
    bool reset_origin;
    bool init_pos;
    bool init_vel;
    bool double_well;
    const char* init_pos_file;
    const char* init_vel_file;
};
void System::init(int amount_in, int timesteps_in, int refresh_in, int threads_in, int dump_in, int reset_steps_in, double dt_in, double temp_in, double mass_in, double friction_in, double density_in, double cutoff_radius_in, bool periodic_in, bool reset_origin_in, bool init_pos_in, bool init_vel_in, const char* init_pos_file_in, const char* init_vel_file_in, bool double_well_in, double difference_in, double second_well_in, double height_in){
// System initializer
  amount = amount_in;
  if (amount == 0) amount = 512;
  timesteps = timesteps_in;
  if (timesteps == 0) timesteps = 1;
  refresh = refresh_in;
  if (refresh == 0) refresh = 1;
  dump = dump_in;
  if (dump == 0) dump = 1;
  reset_steps = reset_steps_in;
  if (reset_steps == 0) reset_steps = 1;
  threads = threads_in;
  if (threads == 0) threads = 1;
  dt = dt_in; 
  if (dt == 0) dt = 0.001;
  init_temp = temp_in;
  if (init_temp == 0) init_temp = 1.;
  current_temp = init_temp;
  mass = mass_in;
  if (mass == 0) mass = 1;
  friction = friction_in;
  density = density_in;
  if (density ==0) density = 0.6;
  region_length = pow((amount*mass)/density,1./3.);
  volume = region_length*region_length*region_length;
  if (region_length == 0) region_length = 100;
  cutoff_radius = cutoff_radius_in;
  if (cutoff_radius == 0) cutoff_radius = region_length/2.;
  periodic = periodic_in;
  reset_origin = reset_origin_in;
  init_pos = init_pos_in;
  init_vel = init_vel_in;
  init_pos_file = init_pos_file_in;
  init_vel_file = init_vel_file_in;
  double_well = double_well_in;
  difference = difference_in;
  second_well = second_well_in;
  height = height_in;
  particles.init(amount, threads, refresh, reset_steps, dt, init_temp, mass, friction, region_length, cutoff_radius, difference, second_well, height, periodic, reset_origin, init_pos, init_vel, double_well, init_pos_file, init_vel_file);
}
double System::kinetic_energy(void){
// Returns the kinetic energy of the system
  double kinetic=0;
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for reduction(+:kinetic)
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double vel = particles.velocities[3*i+j];
      kinetic += vel*vel;
    }
  }
  #pragma omp barrier
  return (mass*kinetic)/2.;
}
double System::potential_energy(void){
// Returns the potential energy of the system
  return particles.potential_energy;
}
double System::temperature(void){
// Returns the current temperature which the system *should* have
  return current_temp;
}
double System::actual_temperature(void){
// Returns the *actual* temperature of the system
  return (kinetic_energy()*2.)/(3.*amount);
}
double System::pressure(void){
// Returns the pressure
  return (1./(3*volume))*(2*kinetic_energy() + particles.virial);
}
void System::update_neighborlist(void){
// Updates the neighborlist
  particles.update_neighborlist();
}
void System::update_forces(void){
// Updates the forces
  particles.update_forces();
}
void System::update_positions(void){
// Updates the positions
  if (friction!=0){
    particles.update_positions();
  }
  else{
    particles.update_positions_verlet();
  }
}
void System::update_velocities(void){
// Updates the velocities
  if (friction!=0){
    particles.update_velocities();
  }
  else{
    particles.update_velocities_verlet();
  }
}
void System::update_rnd(void){
// Recalculates the normally distributed values for the Langevin integrator
  if (friction!=0){
    particles.calc_rnd();
  }
}
void System::set_temperature(double temp_in){
// Sets new temperature and velocities for the system
  current_temp = temp_in;
  particles.sigma = sqrt((2*temp_in*friction)/mass);
  particles.current_temp = temp_in;
}

