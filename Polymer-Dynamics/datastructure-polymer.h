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

double distance_non_periodic(double * xi, double * xj, double * dr){
// Returns the distance between two points in three dimensions
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
    double * positions;
    double * velocities;
    double * forces;
    double * rnd;
    int * boxes;
    int ** neighborlist;
    void init(int amount_in, int threads_in, int refresh_in, double region_length_in, double dt_in, double mass_in, double temp_in, double cutoff_radius_in, double g_in, double distance_in, double epsilon_in, double c_in);
    void init_pos(bool file,const char * filename);
    void init_vel(bool file,const char * filename);
    void update_neighborlist(void);
    void update_forces(void);
    void update_positions_velocity_verlet(void);
    void update_positions_langevin(void);
    void update_velocities_velocity_verlet(void);
    void update_velocities_langevin(void);
    void calc_rnd(void);
  private:
    int amount;
    int threads;
    int called;
    int refresh;
    double region_length;
    double dt;
    double mass;
    double temp;
    double g;
    double sigma;
    double cutoff_radius;
    double verlet_radius;
    double distance;
    double epsilon;
    double c;
};
void Particles::init(int amount_in, int threads_in, int refresh_in, double region_length_in, double dt_in, double mass_in, double temp_in, double cutoff_radius_in, double g_in, double distance_in, double epsilon_in, double c_in){
// Initializer for Particles class
  amount = amount_in;
  threads = threads_in;
  refresh = refresh_in;
  region_length = region_length_in;
  temp = temp_in;
  dt = dt_in;
  mass = mass_in;
  cutoff_radius = cutoff_radius_in;
  verlet_radius = cutoff_radius + 2*refresh*dt;
  g = g_in;
  distance = distance_in;
  epsilon = epsilon_in;
  c = c_in;
  sigma = sqrt((2*temp*g)/mass);
  positions = new double[amount*3];
  velocities = new double[amount*3];
  forces = new double[amount*3];
  boxes = new int[amount*3];
  neighborlist = new int*[amount-1]; // todo: allocate memory based on expected neighbors to save memory
  for (int i=0; i<amount; ++i){
    neighborlist[i] = new int[amount];
  }
  if (g!=0){
    rnd = new double[amount*6];
  }
  called = 0;
}
void Particles::update_neighborlist(void){
  if (called%refresh==0){
    double dr[3];
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<(amount-1); ++i){
      int l = 0;
      for (int j=(i+1); j<amount; ++j){
        double pos_i[3], pos_j[3];
        for (int k=0; k<3; ++k){
          pos_i[k] = positions[3*i+k];
          pos_j[k] = positions[3*j+k];
        }
        double distance_ij = distance_non_periodic(pos_i,pos_j,dr);
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
// Calculates the force based on the Lennard-Jones-Potential, the harmonic potential due to the springs
// and due to the external potential
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
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<(amount-1); ++i){
    int neighbors = neighborlist[i][0];
    for (int j=1; j<=neighbors; ++j){
      double pos_i[3], pos_j[3], dr[3];
      int index_j = neighborlist[i][j];
      for (int k=0; k<3; ++k){
        pos_i[k] = positions[3*i+k];
        pos_j[k] = positions[3*index_j+k];
      }
      double d2 = distance_non_periodic(pos_i,pos_j,dr);
      double di2 = 1./d2;
      double di6 = di2*di2*di2;
      double force = 48.*di6*(di6-0.5)*di2;
      for (int k=0; k<3; ++k){ 
        #pragma omp atomic   
          forces[3*i+k] += force*dr[k];
        #pragma omp atomic
          forces[3*index_j+k] -= force*dr[k];
      }
    }
  }
  #pragma omp barrier
  //////////////////////////
  /// HARMONIC POTENTIAL ///
  //////////////////////////  
  double k = (3.*temp)/(distance*distance);
  for (int i=0; i<3; ++i){
    forces[i] += k*(positions[i+3] - positions[i]);
    forces[3*amount-3+i] += k*(positions[3*amount-6+i] - positions[3*amount-3+i]);
  }
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=1; i<(amount-1); ++i){
    for (int j=0; j<3; ++j){
      forces[3*i+j] += k*(positions[3*i+3+j] - 2*positions[3*i+j] + positions[3*i-3+j]);
    }
  }
  #pragma omp barrier
  //////////////////////////////////////
  /// EXTERNAL DOUBLE WELL POTENTIAL ///
  //////////////////////////////////////
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    double x = positions[3*i];
    forces[3*i] += -((4.*epsilon*x)/(c*c*c*c))*(x*x - c*c);
  }
}
void Particles::update_positions_velocity_verlet(void){
// Updates the postions according to the velocity verlet integrator
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double pos = positions[3*i+j] + dt*velocities[3*i+j];
      positions[3*i+j] = pos;
    }
  }
  #pragma omp barrier
}
void Particles::update_positions_langevin(void){
// Updates the positions according to the Langevin integrator
  double s3dt = sqrt(dt*dt*dt);
  double si12 = 1./sqrt(12.);
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      double pos = positions[3*i+j] + dt*velocities[3*i+j] + s3dt*sigma*si12*rnd[6*i+j+3];
      positions[3*i+j] = pos;
    }
  }
}
void Particles::update_velocities_velocity_verlet(void){
// Updates the velocity according to the velocity verlet integrator
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
}
void Particles::update_velocities_langevin(void){
// Updates the velocity according to the Langevin integrator
  double sdt = sqrt(dt);
  double s3dt = sqrt(dt*dt*dt);
  double si3 = 1./sqrt(3.);
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    for (int j=0; j<3; ++j){
      velocities[3*i+j] += 0.5*dt*forces[3*i+j]/mass - 0.5*dt*g*velocities[3*i+j]  + 0.5*sdt*sigma*rnd[6*i+j] - 0.125*dt*dt*g*(forces[3*i+j]/mass - g*velocities[3*i+j]) - 0.25*s3dt*g*sigma*(0.5*rnd[6*i+j] + si3*rnd[6*i+j+3]);
    }
  }
  #pragma omp barrier
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
void Particles::init_pos(bool file, const char * filename){
// Either reads a file containing the initial positions or sets them in a row equally spaced
  if (!file){
    double sig = distance/sqrt(3);
    positions[0] = -distance;
    positions[1] = 0;
    positions[2] = 0;
    for (int i=1; i<amount; ++i){
      for (int j=0; j<3; ++j){
        double u = (double)rand()/RAND_MAX;
        double v = (double)rand()/RAND_MAX;
        positions[3*i+j] = positions[3*i+j-3] + sig*sqrt(-2*log(u))*cos(2*PI*v);
      }
    }
  }
  else{
    ifstream init_pos_file (filename);
    if (init_pos_file.is_open()){
      // Reads from a .xyz file
      // first row: amount of particles
      // second row: description
      // following rows: particle positions
      int amount_file;
      int i=0;
      string description;
      string atom;
      double x,y,z;
      init_pos_file >> amount_file;
      init_pos_file >> description;
      if (amount == amount_file){
        while (init_pos_file >> atom >> x >> y >> z){
          for (int j=0; j<amount; ++j){
            boxes[3*i+j] = 0;
          }
          positions[3*i] = x;
          positions[3*i+1] = y;
          positions[3*i+2] = z;
        }
      }
      else{
        std::cerr << "The amount of particles inside the file does not match the specified amount by the command line" << std::endl;
        exit(EXIT_FAILURE);
      }
      init_pos_file.close();
    }
    else{
      perror(filename);
    } 
  }
}
void Particles::init_vel(bool file, const char * filename){
// Either reads a file containing the initial velocities or draw them from a Boltzmann-distribution
// and removes net momentum to avoid a center of mass drift
  double mom_x=0, mom_y=0, mom_z=0,u,v,r;
  if (!file){
    for (int i=0; i<amount; ++i){
      for (int j=0; j<3; ++j){
        u = (double)rand()/RAND_MAX;
        v = (double)rand()/RAND_MAX;
        r = sqrt(temp/mass)*sqrt(-2*log(u))*cos(2*PI*v);
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
    ifstream init_vel_file (filename);
    if (init_vel_file.is_open()){
      // Reads from a .xyz file
      // first row: amount of particles
      // second row: description
      // following rows: particle velocities
      int amount_file;
      int i=0;
      string description;
      string atom;
      double x,y,z;
      init_vel_file >> amount_file;
      init_vel_file >> description;
      if (amount == amount_file){
        while (init_vel_file >> atom >> x >> y >> z){
          velocities[3*i] = x;
          velocities[3*i+1] = y;
          velocities[3*i+2] = z;
        }
      }
      else{
        std::cerr << "The amount of particles inside the file does not match the specified amount by the command line" << std::endl;
        exit(EXIT_FAILURE);
      }
      init_vel_file.close();
    }
    else{
      perror(filename);
    } 
  }
}


class System{
  public:
    System(int amount_in, int timesteps_in, int refresh_in, int threads_in, double dt_in, double temp_in, double mass_in, double g_in, double distance_in, double region_length_in, double cutoff_radius_in, double seed_in, double epsilon_in, double c_in, bool file_pos_out, const char * filename_pos_out, bool file_pos_in, const char * filename_pos_in, bool file_vel_in, const char * filename_vel_in);
    void update_neighborlist(void);
    void update_forces(void);
    void update_velocities(void);
    void update_positions(void);
    void update_rnd(void);
    double potential_energy(void);
    double kinetic_energy(void);
    double temperature(void);
    double pressure(void);
    void print_init(void);
    void print_forces(void);
    void print_neighborlist(void);
    void print_positions(void);
    void print_velocities(void);
    void print_total_momentum(void);
    void print_m_to_n(int m, int n);
  private:
    Particles particles;
    int amount;
    int refresh;
    int threads;
    int timesteps;
    double total_time;
    double dt;
    double temp;
    double mass;
    double distance;
    double region_length;
    double cutoff_radius;
    double g;
    double seed;
    double epsilon;
    double c;
};
System::System(int amount_in, int timesteps_in, int refresh_in, int threads_in, double dt_in, double temp_in, double mass_in, double g_in, double distance_in, double region_length_in, double cutoff_radius_in, double seed_in, double epsilon_in, double c_in, bool file_pos_out, const char * filename_pos_out, bool file_pos_in, const char * filename_pos_in, bool file_vel_in, const char * filename_vel_in){
// System initializer
  amount = amount_in;
  refresh = refresh_in;
  timesteps = timesteps_in;
  threads = threads_in;
  dt = dt_in;
  temp = temp_in;
  mass = mass_in;
  g = g_in;
  distance = distance_in;
  seed = seed_in;
  epsilon = epsilon_in;
  c = c_in;
  srand(seed);
  region_length = region_length_in;
  cutoff_radius = cutoff_radius_in;
  total_time = 0;
  particles.init(amount,threads,refresh,region_length,dt,mass,temp,cutoff_radius,g,distance,epsilon,c);
  particles.init_pos(file_pos_in,filename_pos_in);
  particles.init_vel(file_vel_in,filename_vel_in);
}
void System::update_neighborlist(void){
// Updates the neighborlist
  particles.update_neighborlist();
}
void System::update_forces(void){
// Updates the forces
  particles.update_forces();
}
void System::update_velocities(void){
// Updates the velocities
  if (g!=0){
    particles.update_velocities_langevin();
  }
  else{
    particles.update_velocities_velocity_verlet();
  }
}
void System::update_positions(void){
// Updates the positions
  total_time += dt;
  if (g!=0){
    particles.update_positions_langevin();
  }
  else{
    particles.update_positions_velocity_verlet();
  }
}
void System::update_rnd(void){
// Recalculates the normally distributed values for the Langevin integrator
  particles.calc_rnd();
}
double System::kinetic_energy(void){
// Returns the kinetic energy
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
double System::temperature(void){
// Returns the temperature
  return (kinetic_energy()*2.)/(3.*amount);
}
void System::print_init(void){
// Prints out initial data	
  std::cerr << "# particles:\t\t" << amount << std::endl;
  std::cerr << "refresh rate:\t\t" << refresh << std::endl;
  std::cerr << "temperature:\t\t"<< temp << std::endl;
  std::cerr << "box sidelength:\t\t"<< region_length << std::endl;
  std::cerr << "cutoff radius:\t\t"<< cutoff_radius << std::endl;
  std::cerr << "mass:\t\t\t" << mass << std::endl;
  std::cerr << "friction:\t\t" << g << std::endl;
  std::cerr << "# timesteps:\t\t" << timesteps << std::endl;
  std::cerr << "timestep:\t\t" << dt << std::endl;
  std::cerr << "mean bond length:\t" << distance << std::endl;
  std::cerr << "potential depth:\t" << epsilon << std::endl;
  std::cerr << "potential root:\t\t" << c << std::endl;
  std::cerr << "# threads:\t\t" << threads << "\n" << std::endl;
}
void System::print_forces(void){
// Prints the forces
  printf("t = %4.4f\n",total_time);
  for (int i=0; i<amount; ++i){
    printf("force (%i): %4.4f\t%4.4f\t%4.4f\n",i,particles.forces[3*i],particles.forces[3*i+1],particles.forces[3*i+2]);
  }
  printf("\n");
}
void System::print_neighborlist(void){
// Prints the neighborlist
  printf("t = %4.4f\n",total_time);
  for (int i=0; i<amount; ++i){
    int neighbors  = particles.neighborlist[i][0];
    printf("neighbors (%i): ", i);
    for (int j=1; j<=neighbors; ++j){
      printf("%i ",particles.neighborlist[i][j]);
    }
    printf("\n");
  }
}
void System::print_positions(void){
// Prints the positions (xyz compatible format)
  printf("%i\n",amount);
  printf("Particle positions\n");
  for (int i=0; i<amount; ++i){
    printf("Ar\t%4.5f\t%4.5f\t%4.5f\n",particles.positions[3*i],particles.positions[3*i+1],particles.positions[3*i+2]);
  }
}
void System::print_velocities(void){
// Prints the velocities
  printf("t = %4.4f\n",total_time);
  for (int i=0; i<amount; ++i){
    printf("velocity (%i): %4.4f\t%4.4f\t%4.4f\n",i,particles.velocities[3*i],particles.velocities[3*i+1],particles.velocities[3*i+2]); 
  }
  printf("\n");
}
void System::print_total_momentum(void){
// Prints the total momentum
  double mom_x=0, mom_y=0, mom_z=0;
  for (int i=0; i<amount; ++i){
    mom_x += mass*particles.velocities[3*i];
    mom_y += mass*particles.velocities[3*i+1];
    mom_z += mass*particles.velocities[3*i+2];
  }
  printf("%4.4f\t%4.4f\t%4.4f\n",mom_x,mom_y,mom_z); 
}
void System::print_m_to_n(int m, int n){
// Prints the end to end distance between the m-th and the n-th bead
  double m_to_n[3];
  for (int i=0; i<3; ++i){
    m_to_n[i] = 0;
  }
  for (int i=m; i<n; ++i){
    for (int j=0; j<3; ++j){
      m_to_n[j] += particles.positions[3*i+j] - particles.positions[3*i-3+j];
    }
  }
  for (int i=0; i<3; ++i){
    printf("%4.4f\t", m_to_n[i]);
  }
  printf("\n");
}
