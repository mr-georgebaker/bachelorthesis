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
// for periodic boundary conditions
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

class Particles{
  public:
    double * positions;
    double * velocities;
    double * forces;
    double * rnd;
    int * boxes;
    void init(int amount_in, int threads_in, double region_length_in, double dt_in, double mass_in, double temp_in, double cutoff_radius_in, double g_in, double distance_in, double epsilon_in, double c_in, bool excluded_in);
    void init_pos(bool file,const char * filename);
    void init_vel(bool file,const char * filename);
    void update_forces(void);
    void update_positions_langevin(void);
    void update_velocities_langevin(void);
    void calc_rnd(void);
  private:
    int amount;
    int threads;
    int called;
    double region_length;
    double dt;
    double mass;
    double temp;
    double g;
    double sigma;
    double cutoff_radius;
    double rc2;
    double distance;
    double epsilon;
    double c;
    double spring;
    bool excluded;
    double wrap_positions(double pos_in, int index);
};
void Particles::init(int amount_in, int threads_in, double region_length_in, double dt_in, double mass_in, double temp_in, double cutoff_radius_in, double g_in, double distance_in, double epsilon_in, double c_in, bool excluded_in){
// Initializer for Particles class
  amount = amount_in;
  threads = threads_in;
  temp = temp_in;
  dt = dt_in;
  mass = mass_in;
  cutoff_radius = cutoff_radius_in;
  rc2 = cutoff_radius*cutoff_radius;
  g = g_in;
  distance = distance_in;
  region_length = region_length_in;
  epsilon = epsilon_in;
  c = c_in;
  sigma = sqrt((2*temp*g)/mass);
  positions = new double[amount*3];
  velocities = new double[amount*3];
  forces = new double[amount*3];
  boxes = new int[amount*3];
  spring = (3.*temp)/(distance*distance);
  if (g!=0){
    rnd = new double[amount*6];
  }
  excluded = excluded_in;
}
void Particles::update_forces(void){
// Calculates the forces based on the Lennard-Jones-Potential, the harmonic potential due to the springs
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
  if (excluded){
    #ifdef _OPENMP
      omp_set_num_threads(threads);
      #pragma omp parallel for
    #endif
    for (int i=0; i<(amount-1); ++i){
      for (int j=i+1; j<amount; ++j){
        double pos_i[3], pos_j[3], dr[3];
        for (int k=0; k<3; ++k){
          pos_i[k] = positions[3*i+k];
          pos_j[k] = positions[3*j+k];
        }
        double distance_ij = distance_periodic(pos_i,pos_j,dr,region_length);
        if (distance_ij < rc2){
          double di2 = 1./distance_ij;
          double di6 = di2*di2*di2;
          double force = 48.*di6*(di6-0.5)*di2;
          for (int k=0; k<3; ++k){
            #pragma omp atomic   
              forces[3*i+k] += force*dr[k];
            #pragma omp atomic
              forces[3*j+k] -= force*dr[k];
          }
        }
      }
    }  
  }
  #pragma omp barrier
  //////////////////////////
  /// HARMONIC POTENTIAL ///
  //////////////////////////  
  for (int i=0; i<3; ++i){
    forces[i] += spring*(positions[i+3] - positions[i]);
    forces[3*amount-3+i] += spring*(positions[3*amount-6+i] - positions[3*amount-3+i]);
  }
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=1; i<(amount-1); ++i){
    for (int j=0; j<3; ++j){
      forces[3*i+j] += spring*(positions[3*i+3+j] - 2*positions[3*i+j] + positions[3*i-3+j]);
    }
  }
  #pragma omp barrier
  //////////////////////////////////////
  /// EXTERNAL DOUBLE WELL POTENTIAL ///
  //////////////////////////////////////
  /*
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for
  #endif
  for (int i=0; i<amount; ++i){
    double x = positions[3*i];
    forces[3*i] += -((4.*epsilon*x)/(c*c*c*c))*(x*x - c*c);
  }
  */
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
      positions[3*i+j] = wrap_positions(pos,3*i+j);
    }
  }
}
double Particles::wrap_positions(double pos_in, int index){
// Wraps the position to keep the particle inside the box
  if (pos_in > region_length){
    pos_in -= region_length;
    boxes[index] += 1; 
  }
  else if (pos_in < 0){
    pos_in += region_length;
    boxes[index] -= 1;
  }
  return pos_in;
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
// Either reads a file containing the initial positions or sets them spaced accoring to a Gaussian distribution
  if (!file){
    double r;
    double sig = distance/sqrt(3);
    positions[0] = region_length/2-(amount*distance)/2.;
    positions[1] = region_length/2.;
    positions[2] = region_length/2.;
    for (int i=1; i<amount; ++i){
      for (int j=0; j<3; ++j){
        do{
          double u = (double)rand()/RAND_MAX;
          double v = (double)rand()/RAND_MAX;
          r = sig*sqrt(-2*log(u))*cos(2*PI*v);
        }while(r<=0.9);
        positions[3*i+j] = positions[3*i+j-3] + r;
        boxes[3*i+j] = 0;
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
          ++i;
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
    System(int amount_in, int timesteps_in, int threads_in, double dt_in, double temp_in, double mass_in, double g_in, double distance_in, double region_length_in, double cutoff_radius_in, double seed_in, double epsilon_in, double c_in, bool file_pos_out, const char * filename_pos_out, bool file_pos_in, const char * filename_pos_in, bool file_vel_in, const char * filename_vel_in, bool excluded_in);
    void update_forces(void);
    void update_velocities(void);
    void update_positions(void);
    void update_rnd(void);
    void print_init(int dump);
    void print_positions(void);
    void print_m_to_n(int m, int n);
  private:
    Particles particles;
    int amount;
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
    bool excluded;
};
System::System(int amount_in, int timesteps_in, int threads_in, double dt_in, double temp_in, double mass_in, double g_in, double distance_in, double region_length_in, double cutoff_radius_in, double seed_in, double epsilon_in, double c_in, bool file_pos_out, const char * filename_pos_out, bool file_pos_in, const char * filename_pos_in, bool file_vel_in, const char * filename_vel_in, bool excluded_in){
// System initializer
  amount = amount_in;
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
  excluded = excluded_in;
  srand(seed);
  region_length = region_length_in;
  if (region_length==0) region_length = 50*amount*distance;
  cutoff_radius = cutoff_radius_in;
  if (cutoff_radius==0) cutoff_radius = region_length/2.;
  total_time = 0;
  particles.init(amount,threads,region_length,dt,mass,temp,cutoff_radius,g,distance,epsilon,c,excluded);
  particles.init_pos(file_pos_in,filename_pos_in);
  particles.init_vel(file_vel_in,filename_vel_in);
}
void System::update_forces(void){
// Updates the forces
  particles.update_forces();
}
void System::update_velocities(void){
// Updates the velocities
  particles.update_velocities_langevin();
}
void System::update_positions(void){
// Updates the positions
  total_time += dt;
  particles.update_positions_langevin();
}
void System::update_rnd(void){
// Recalculates the normally distributed values for the Langevin integrator
  particles.calc_rnd();
}
void System::print_init(int dump){
// Prints out initial data	
  std::cerr << "# particles:\t\t" << amount << std::endl;
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
  std::cerr << "# threads:\t\t" << threads << std::endl;
  std::cerr << "excluded volume:\t" << excluded << std::endl;
  std::cerr << "# steps before dumping values:\t" << dump << std::endl;
  std::cerr << std::endl;
  std::cerr << "relaxation time of first mode:\t\t" << (g*amount*amount*distance*distance)/(3*PI*PI) << std::endl;
  std::cerr << "steps needed to achieve relaxation:\t" << ceil(((g*amount*amount*distance*distance)/(3*PI*PI))*(1/dt)) << "\n" << std::endl;
}
void System::print_positions(void){
// Prints the positions (.xyz compatible format)
  printf("%i\n",amount);
  printf("Particle positions\n");
  for (int i=0; i<amount; ++i){
    printf("Ar\t%4.5f\t%4.5f\t%4.5f\n",particles.positions[3*i],particles.positions[3*i+1],particles.positions[3*i+2]);
  }
}
void System::print_m_to_n(int m, int n){
// Prints the squared distance between the m-th and n-th bead
  double pos_n[3], pos_m[3], dr[3];
  for (int i=0; i<3; ++i){
    pos_n[i] = particles.positions[3*n-3+i];
    pos_m[i] = particles.positions[3*m-3+i];
  }
  double distance_mn = distance_periodic(pos_n,pos_m,dr,region_length);
  printf("%4.4f\n", distance_mn);
}
