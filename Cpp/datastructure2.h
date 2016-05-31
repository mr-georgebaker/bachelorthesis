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
    double potential_energy;
    double virial;
    double * positions;
    double * velocities;
    double * forces;
    double * rnd;
    int * boxes;
    int ** neighborlist;
    void init(int amount_in,int threads_in,int refresh_in,double region_length_in,double dt_in,double mass_in,double temp_in,double cutoff_radius_in,double g_in);
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
    double wrap_positions(double pos_in);
};
void Particles::init(int amount_in,int threads_in,int refresh_in,double region_length_in,double dt_in,double mass_in,double temp_in,double cutoff_radius_in,double g_in){
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
        double distance_ij = distance_periodic(pos_i,pos_j,dr,region_length);
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
// Calculates the force and potential energy based on the Lennard-Jones potential 
  double pot_energy = 0, vir = 0; 
  double dvi2 = 1./(verlet_radius*verlet_radius);
  double dvi6 = dvi2*dvi2*dvi2;
  double potential_energy_rc = 4.*dvi6*(dvi6-1.);
  //double force_rc = 48.*dvi6*(dvi6-0.5)*dvi2;
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
  #ifdef _OPENMP
    omp_set_num_threads(threads);
    #pragma omp parallel for reduction(+:vir,pot_energy)
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
      double d2 = distance_periodic(pos_i,pos_j,dr,region_length);
      double di2 = 1./d2;
      double di6 = di2*di2*di2;
      double force = 48.*di6*(di6-0.5)*di2;//-force_rc; 
      vir += force;
      for (int k=0; k<3; ++k){ 
        #pragma omp atomic   
          forces[3*i+k] += force*dr[k];
        #pragma omp atomic
          forces[3*index_j+k] -= force*dr[k];
      }
      pot_energy += 4.*di6*(di6-1.) - potential_energy_rc;// + force_rc*(verlet_radius - sqrt(d2));
    }
  }
  #pragma omp barrier
  virial = vir;
  potential_energy = pot_energy;
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
      positions[3*i+j] = wrap_positions(pos);
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
      positions[3*i+j] = wrap_positions(pos);
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
// Either reads a file containing the initial positions or sets them on a crystal like grid
  int ix=0, iy=0, iz=0, spacing=1;
  if (!file){
    while ((spacing*spacing*spacing)<amount) spacing++;
    for (int i=0; i<amount; ++i){
      positions[3*i] = ((ix+0.5)*region_length)/spacing;
      positions[3*i+1] = ((iy+0.5)*region_length)/spacing;
      positions[3*i+2] = ((iz+0.5)*region_length)/spacing;
      for (int j=0; j<3; ++j){
        boxes[3*i+j] = 0;
      }
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
    System(int amount_in, int refresh_in, int timesteps_in, int threads_in, double dt_in, double temp_in, double mass_in, double rho_in, double radius_in, double g_in, double seed_in, bool init_pos, const char * init_pos_file, bool init_vel, const char * init_vel_file);
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
    void print_energy(void);
    void print_kinetic_energy(void);
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
    double rho;
    double region_length;
    double cutoff_radius;
    double g;
    double seed;
    double avk;
    double avk2;
};
System::System(int amount_in, int refresh_in, int timesteps_in, int threads_in, double dt_in, double temp_in, double mass_in, double rho_in, double cutoff_radius_in, double g_in, double seed_in, bool init_pos, const char * init_pos_file, bool init_vel, const char * init_vel_file){
// System initializer
  amount = amount_in;
  refresh = refresh_in;
  timesteps = timesteps_in;
  threads = threads_in;
  dt = dt_in;
  temp = temp_in;
  mass = mass_in;
  rho = rho_in;
  g = g_in;
  seed = seed_in;
  srand(seed);
  region_length = pow((amount*mass)/rho,1./3.);
  cutoff_radius = cutoff_radius_in;
  if(cutoff_radius==0) cutoff_radius = region_length/2.;
  avk = 0;
  avk2 = 0;
  total_time = 0;
  particles.init(amount,threads,refresh,region_length,dt,mass,temp,cutoff_radius,g);
  particles.init_pos(init_pos,init_pos_file);
  particles.init_vel(init_vel,init_vel_file);
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
double System::potential_energy(void){
// Returns the potential energy
  return particles.potential_energy;
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
double System::pressure(void){
// Returns the pressure
  return (1./(3.*region_length*region_length*region_length))*((2./3.)*kinetic_energy() + particles.virial);
}
void System::print_init(void){
// Prints out initial data	
  std::cerr << "# particles:\t" << amount << std::endl;
  std::cerr << "refresh rate:\t" << refresh << std::endl;
  std::cerr << "temperature:\t"<< temp << std::endl;
  std::cerr << "box sidelength:\t"<< region_length << std::endl;
  std::cerr << "cutoff radius:\t"<< cutoff_radius << std::endl;
  std::cerr << "mass:\t\t" << mass << std::endl;
  std::cerr << "density:\t" << rho << std::endl;
  std::cerr << "friction:\t" << g << std::endl;
  std::cerr << "# timesteps:\t" << timesteps << std::endl;
  std::cerr << "timestep:\t" << dt << std::endl;
  std::cerr << "# threads:\t" << threads << "\n" << std::endl;
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
    printf("Ar\t%4.4f\t%4.4f\t%4.4f\n",particles.positions[3*i],particles.positions[3*i+1],particles.positions[3*i+2]);
  }
}
void System::print_velocities(void){
// Prints the velocities
  for (int i=0; i<amount; ++i){
    printf("%4.4f\t%4.4f\t%4.4f\n",particles.velocities[3*i],particles.velocities[3*i+1],particles.velocities[3*i+2]); 
  }
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
void System::print_energy(void){
// Prints the kinetic, potential and total energy
  double kin = kinetic_energy();
  double pot = potential_energy();
  printf("%4.4f\t%4.2f\t%4.2f\t%4.2f\n",total_time,kin,pot,kin+pot);
}
void System::print_kinetic_energy(void){
// Prints the kinetic energy and its standard deviation
  double kin = kinetic_energy();
  int called = total_time/dt;
  avk += kin;
  avk2 += kin*kin;
  double var = (avk2 - (avk*avk/called))/called;
  printf("%4.2f\t%4.2f\n",kin,sqrt(var));
}
