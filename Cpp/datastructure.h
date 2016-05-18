#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

#define PI 3.14159265

using namespace std;

class Vec{
  public:
    double vec[3];
};

double distance_periodic(Vec xi, Vec xj, double region, Vec * dr){
// Returns the squared distance between two points in three dimensions
// for periodic boundary conditions
  double d=0, r=region/2., c;
  c = xi.vec[0] - xj.vec[0];
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
  dr->vec[0] = c ;
  c = xi.vec[1] - xj.vec[1];
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
  dr->vec[1] = c ;
  c = xi.vec[2] - xj.vec[2];
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
  dr->vec[2] = c ;
  return d;
}

class Atom{
  public:
    int index;
    double mass;
    Vec box;
    Vec position;
    Vec velocity;
    Vec force;
};

class XiEta{
  public:
    Vec ** xieta;
    void init(int amount_in, int num_threads_in);
    void calc(void);
  private:
    int amount;
    int num_threads;
};
void XiEta::init(int amount_in, int num_threads_in){
// Initializer - Allocates memory and sets initial variables
  amount = amount_in;
  num_threads = num_threads_in;
  xieta = new Vec*[amount];
  for (int i=0; i<amount; ++i){
    xieta[i] = new Vec[2];
  }
}
void XiEta::calc(void){
// Calculates normally distributed values with mu=0 and sigma=1 
// based on the box muller method
  double u,v;
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for private(u,v)
  #endif
    for (int i=0; i<amount; i++){
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].vec[0] = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].vec[0] = sqrt(-2*log(u))*sin(2*PI*v);
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].vec[1] = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].vec[1] = sqrt(-2*log(u))*sin(2*PI*v);
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].vec[2] = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].vec[2] = sqrt(-2*log(u))*sin(2*PI*v);
    }
}

class Particles{
  public:
    Atom * particlelist;
    double potential_energy;
    void init(int amount_in, int num_threads_in, double region_in, double mass_in, double temp_in, double dt_in, double radius_in);
    void init_positions(bool file, const char * filename);
    void init_velocities(bool file, const char * filename);
    void update_force_lennard_jones(int ** neighborlist);
    void update_velocity_verlet(void);
    void update_position_verlet(void);
  private:
    int amount;
    int num_threads;
    double region; 
    double mass;
    double temp;
    double dt;
    double radius;
    double wrap_force(double d_in);
    double wrap_position(double coord_in);
};
void Particles::init(int amount_in, int num_threads_in, double region_in, double mass_in, double temp_in, double dt_in, double radius_in){
// Initializes the values for the Particles class and allocates memory for the particlelist
  amount = amount_in;
  num_threads = num_threads_in;
  region = region_in;
  mass = mass_in;
  temp = temp_in;
  dt = dt_in;
  radius = radius_in;
  particlelist = new Atom[amount];
}
void Particles::init_positions(bool file=false, const char * filename="inital_positions.xyz"){
// Either reads a file containing the initial positions or sets them on a crystal like grid
  int ix=0, iy=0, iz=0, spacing=1;
  if (!file){
    while ((spacing*spacing*spacing)<amount) spacing++;
    for (int i=0; i<amount; ++i){
      particlelist[i].mass = mass;
      particlelist[i].index = i;
      particlelist[i].box.vec[0] = 0;
      particlelist[i].box.vec[1] = 0;
      particlelist[i].box.vec[2] = 0;
      particlelist[i].position.vec[0] = ((ix+0.5)*region)/spacing;
      particlelist[i].position.vec[1] = ((iy+0.5)*region)/spacing;
      particlelist[i].position.vec[2] = ((iz+0.5)*region)/spacing;
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
    FILE *initial_pos = fopen(filename, "r");
    if (initial_pos != NULL){
      // Todo: Read positions from a file
      fclose(initial_pos);
    }
    else{
      perror(filename);
    }
  }
}
void Particles::init_velocities(bool file=false, const char * filename="inital_velocities.xyz"){
// Either reads a file containing the initial velocities or draw them from a Boltzmann-distribution
// and removes net momentum to avoid a center of mass drift
  double mom_x, mom_y, mom_z;
  double u, v, r;
  if (!file){
    #ifdef _OPENMP
      omp_set_num_threads(num_threads);
      #pragma omp parallel for private(u,v,r) reduction(+:mom_x,mom_y,mom_z) 
    #endif
      for (int i=0; i<amount; ++i){
        u = (double)rand()/RAND_MAX;
        v = (double)rand()/RAND_MAX;
        r = sqrt(temp/mass)*sqrt(-2*log(u))*cos(2*PI*v);
        particlelist[i].velocity.vec[0] = r;
        mom_x += r;
        u = (double)rand()/RAND_MAX;
        v = (double)rand()/RAND_MAX;
        r = sqrt(temp/mass)*sqrt(-2*log(u))*cos(2*PI*v);
        particlelist[i].velocity.vec[1] = r;
        mom_y += r;
        u = (double)rand()/RAND_MAX;
        v = (double)rand()/RAND_MAX;
        r = sqrt(temp/mass)*sqrt(-2*log(u))*cos(2*PI*v);
        particlelist[i].velocity.vec[2] = r;
        mom_z += r;
      }
    #pragma omp barrier
    mom_x /= amount;
    mom_y /= amount;
    mom_z /= amount;
    #ifdef _OPENMP
      omp_set_num_threads(num_threads);
      #pragma omp parallel for
    #endif
      for (int i=0; i<amount; ++i){
        particlelist[i].velocity.vec[0] -= mom_x;
        particlelist[i].velocity.vec[1] -= mom_y;
        particlelist[i].velocity.vec[2] -= mom_z;
      } 
    #pragma omp barrier
  }
  else{
    FILE *initial_pos = fopen(filename, "r");
    if (initial_pos != NULL){
      // Todo: Read velocities from a file
      fclose(initial_pos);
    }
    else{
      perror(filename);
    }
  }
}
void Particles::update_force_lennard_jones(int ** neighborlist){
// Updates the force based on the neighborlist and calculates
// the potential energy for the lennard jones potential
  double pot_energy = 0;
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; ++i){
      particlelist[i].force.vec[0] = 0;
      particlelist[i].force.vec[1] = 0;
      particlelist[i].force.vec[2] = 0;
    }
  #pragma omp barrier
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for reduction(+:pot_energy)
  #endif
    for (int i=0; i<(amount-1); ++i){
      int neighbors = neighborlist[i][0];
      for (int j=1; j<=neighbors; ++j){
        Vec dr;
        int index = neighborlist[i][j];
        Vec pos_i = particlelist[i].position;
        Vec pos_j = particlelist[index].position;
        double d2 = distance_periodic(pos_i,pos_j,region, &dr);
        double di2 = 1./d2;
        double di6 = di2*di2*di2;
        double force = 48.*di6*(di6-0.5)*di2;
        double force_x = force*dr.vec[0];
        double force_y = force*dr.vec[1];
        double force_z = force*dr.vec[2];
        #pragma omp atomic
          particlelist[i].force.vec[0] += force_x;
        #pragma omp atomic
          particlelist[i].force.vec[1] += force_y;
        #pragma omp atomic
          particlelist[i].force.vec[2] += force_z;
        #pragma omp atomic
          particlelist[index].force.vec[0] -= force_x;
        #pragma omp atomic
          particlelist[index].force.vec[1] -= force_y;
        #pragma omp atomic
          particlelist[index].force.vec[2] -= force_z;
        pot_energy += di6*(di6-1.);
      }
    }
  #pragma omp barrier
  potential_energy = 4*pot_energy;
}
double Particles::wrap_force(double d_in){
// Applies periodic boundary conditions to the force
  double r = region/2.;
  if (d_in>r){
    d_in -= region;
  }
  else if (d_in<-r){
    d_in += region;
  }
  return d_in;
}
void Particles::update_velocity_verlet(void){
// Updates the velocity based on the velocity verlet integrator
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; ++i){
      particlelist[i].velocity.vec[0] += 0.5*dt*particlelist[i].force.vec[0]/mass;
      particlelist[i].velocity.vec[1] += 0.5*dt*particlelist[i].force.vec[1]/mass;
      particlelist[i].velocity.vec[2] += 0.5*dt*particlelist[i].force.vec[2]/mass;
    }
  #pragma omp barrier
}
void Particles::update_position_verlet(void){
// Updates the position based on the velocity verlet integrator
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; i++){
      double pos_x = particlelist[i].position.vec[0] + dt*particlelist[i].velocity.vec[0];
      double pos_y = particlelist[i].position.vec[1] + dt*particlelist[i].velocity.vec[1];
      double pos_z = particlelist[i].position.vec[2] + dt*particlelist[i].velocity.vec[2];
      pos_x = wrap_position(pos_x);
      pos_y = wrap_position(pos_y);
      pos_z = wrap_position(pos_z);
      particlelist[i].position.vec[0] = pos_x;
      particlelist[i].position.vec[1] = pos_y;
      particlelist[i].position.vec[2] = pos_z;
    }
}
double Particles::wrap_position(double coord_in){
// Wraps the position to keep the particle inside the box
  if (coord_in > region){
    coord_in -= region;
  }
  else if (coord_in < 0){
    coord_in += region;
  }
  return coord_in;
}

class Neighborlist{
  public:
    int ** neighborlist;
    void init(int amount_in, int refresh_in, int num_threads_in, double dt_in, double region_in, double radius_in);
    void update_neighborlist(Atom * particlelist);
    double tot_time;
  private:
    int amount;
    int refresh;
    int num_threads;
    int called;
    int amount_expected_neighbors;
    int * amount_neighbors;
    double dt;
    double region;
    double radius;
    void allocate_memory(int size);
    int max_neighbors(void);
    int expected_neighbors(double radius_in);
    double verlet_radius(Atom * particlelist);
    double max_vel(Atom * particlelist);
};
void Neighborlist::init(int amount_in, int refresh_in, int num_threads_in, double dt_in, double region_in, double radius_in){
// Initializes the values for the Neighborlist class and allocates memory for the neighborlist
  amount = amount_in;
  refresh = refresh_in;
  num_threads = num_threads_in;
  dt = dt_in;
  region = region_in;
  radius = radius_in;
  called = 0;
  amount_neighbors = new int[amount];
  allocate_memory(amount);
}
void Neighborlist::allocate_memory(int size){
// Allocates memory for the neighborlist according to the input
  neighborlist = new int*[amount];
  for (int i=0; i<amount; ++i){
    neighborlist[i] = new int[size];
  }
}
int Neighborlist::expected_neighbors(double radius_in){
// Returns the expected amount of neighbors within the given radius
  double r = radius_in;
  return (((4./3.)*PI*r*r*r)/(region*region*region))*amount;
}
void Neighborlist::update_neighborlist(Atom * particlelist){
// Updates the neighborlist
  if (called%refresh==0){
    double rv = verlet_radius(particlelist) ;
    Vec dr;
    #ifdef _OPENMP
      omp_set_num_threads(num_threads);
      #pragma omp parallel for
    #endif
      for (int i=0; i<(amount-1); ++i){
        int l = 0;
        for (int j=i+1; j<amount; ++j){
          double distance_ij = distance_periodic(particlelist[i].position,particlelist[j].position,region,&dr);
          if (distance_ij<rv*rv){
            neighborlist[i][l+1] = particlelist[j].index;
            l++;          
          }
        }
        amount_neighbors[i] = l;
        neighborlist[i][0] = amount_neighbors[i];
      }
    }
    called++;
    tot_time=called*dt;
  #pragma omp barrier
}
double Neighborlist::verlet_radius(Atom * particlelist){
// Returns the verlet radius (rv = rc + r*v*dt)
  double vmax = max_vel(particlelist);
  return radius + refresh*vmax*dt;
}
double Neighborlist::max_vel(Atom * particlelist){
// Calculates the maximum abolute value of all velocities 
  double c, max1=0, max2=0;  
  for (int i=0; i<amount; ++i){
    c = particlelist[i].velocity.vec[0];
    max1 += c*c;
    c = particlelist[i].velocity.vec[1];
    max1 += c*c;
    c = particlelist[i].velocity.vec[2];
    max1 += c*c;
    if (max1 >= max2){
      max2 = max1;
      max1 = 0;
    }
    else {
      max1 = 0;
    }
  }
  return sqrt(max2);
}
int Neighborlist::max_neighbors(void){
// Returns the maximum amount of neighbors for the particles
  int max1, max2=0;
  for (int i=0; i<amount; ++i){
    max1 = amount_neighbors[i];
    if (max1 > max2){
      max2 = max1;
    }
  }
  return max2;
}

class System{
  public:
    Particles particles;
    System(int amount_in, int refresh_in, int timesteps_in, int num_threads_in, double dt_in, double temp_in, double mass_in, double rho_in, double radius_in, double g_in, double seed_in);
    void update_neighborlist(void);
    void update_force(void);
    void update_velocity(void);
    void update_position(void);
    double kinetic_energy(void);
    void dump_force(void);
    double potential_energy(void);
    void printout(void);
    double time(void);
    Vec total_momentum(void);
  private:
    int amount;
    int refresh;
    int timesteps; 
    int num_threads;
    double dt;
    double temp;
    double mass;
    double rho;
    double region;
    double g;
    double radius;
    double seed;
    double sigma;
    Neighborlist neighborlist;
    XiEta rnd_list;
};
System::System(int amount_in, int refresh_in, int timesteps_in, int num_threads_in, double dt_in, double temp_in, double mass_in, double rho_in, double radius_in, double g_in, double seed_in){
// System initializer
  amount = amount_in;
  refresh = refresh_in;
  timesteps = timesteps_in;
  num_threads = num_threads_in;
  dt = dt_in;
  temp = temp_in;
  mass = mass_in;
  rho = rho_in;
  g = g_in;
  seed = seed_in;
  srand(seed);
  region = pow((amount*mass)/rho,1./3.);
  radius = radius_in;
  if(radius==0) radius = region/2.;
  sigma = sqrt((2*temp*g)/mass);
  particles.init(amount,num_threads,region,mass,temp,dt,radius);
  particles.init_positions();
  particles.init_velocities();
  neighborlist.init(amount,refresh,num_threads,dt,region,radius);
  rnd_list.init(amount,num_threads);
}
void System::update_neighborlist(void){
// Neighborlist update
  neighborlist.update_neighborlist(particles.particlelist);
}
void System::update_force(void){
// Force update
  particles.update_force_lennard_jones(neighborlist.neighborlist);
}
void System::update_velocity(void){
// Velocity update
  particles.update_velocity_verlet();
}
void System::update_position(void){
// Position update
  particles.update_position_verlet();
}
double System::kinetic_energy(void){
// Returns the kinetic energy
  double vel, kinetic=0;
  for (int i=0; i<amount; ++i){
    vel = particles.particlelist[i].velocity.vec[0];
    kinetic += vel*vel;
    vel = particles.particlelist[i].velocity.vec[1];
    kinetic += vel*vel;
    vel = particles.particlelist[i].velocity.vec[2];
    kinetic += vel*vel;
  }
  kinetic = (mass*kinetic)/2.;
  return kinetic;
}

double System::time(void){
// Returns the potential energy
  return neighborlist.tot_time;
}

double System::potential_energy(void){
// Returns the potential energy
  return particles.potential_energy;
}
Vec System::total_momentum(void){
// Returns the total momentum
  Vec momentum;
  momentum.vec[0] = 0;
  momentum.vec[1] = 0;
  momentum.vec[2] = 0;
  for (int i=0; i<amount; ++i){
    momentum.vec[0] += mass*particles.particlelist[i].velocity.vec[0];
    momentum.vec[1] += mass*particles.particlelist[i].velocity.vec[1];
    momentum.vec[2] += mass*particles.particlelist[i].velocity.vec[2];
  }
  return momentum;
}
void System::printout(void){
	
  std::cerr << "# T\t= "<< temp << std::endl;
  std::cerr << "# Box\t= "<< region << std::endl;
  std::cerr << "# Cutoff radius\t= "<< radius << std::endl;
}
void System::dump_force(void){
	for(int i=0;i<amount;i++)
	    std::cerr << "t="<< time() << "force ( "<<i<<" ) : "<<particles.particlelist[i].force.vec[0] << " " <<  particles.particlelist[i].force.vec[1] << " " << particles.particlelist[i].force.vec[2] << "\n" << std::endl;
}
