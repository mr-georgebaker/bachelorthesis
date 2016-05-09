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
    double x;
    double y;
    double z;
};

double distance_periodic(Vec xi, Vec xj, double region){
// Returns the squared distance between two points in three dimensions
// for periodic boundary conditions
  double d=0, r=region/2., c;
  c = xi.x - xj.x;
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
  c = xi.y - xj.y;
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
  c = xi.z - xj.z;
  if (c>r){
    c -= region;
  }
  else if (c<-r){
    c += region;
  }
  d += c*c;
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
    ~XiEta();
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
  double u, v;

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for private(u,v)
  #endif
    for (int i=0; i<amount; i++){
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].x = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].x = sqrt(-2*log(u))*sin(2*PI*v);
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].y = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].y = sqrt(-2*log(u))*sin(2*PI*v);
      u = (double)rand()/RAND_MAX;
      v = (double)rand()/RAND_MAX;
      xieta[i][0].z = sqrt(-2*log(u))*cos(2*PI*v);
      xieta[i][1].z = sqrt(-2*log(u))*sin(2*PI*v);
    }
  #pragma omp barrier
}
XiEta::~XiEta(void){
// Destructor - Deallocates memory
  for (int i=0; i<amount; i++){
    xieta[i] = 0;
    delete[] xieta[i];
  }
  xieta = 0;
  delete[] xieta;
}

class Particles{
  public:
    Atom * particlelist;
    void init(int amount_in, int num_threads_in, double region_in, double m_in);
    void init_positions(void);
    void init_velocities(void);
    void update_force(int ** nlist);
    void update_velocity(Vec ** xieta, double dt, double g, double sigma);
    void update_position(Vec ** xieta, double dt, double sigma);
    ~Particles();
  private:
    int amount;
    int spacing;
    int num_threads;
    double region;
    double mass;
    double wrap_position(double coord_in);
    double rand_range(double min, double max);
};
void Particles::init(int amount_in, int num_threads_in, double region_in, double m_in){
// Sets initial variables and allocates memory
  amount = amount_in;
  region = region_in;
  num_threads = num_threads_in;
  mass = m_in;
  particlelist = new Atom[amount];
}
double Particles::rand_range(double min, double max){
// Returns a uniformally distributed random number between min and max
  double scaled = (double)rand()/RAND_MAX;
  return scaled*(max-min)+min;
}
void Particles::init_positions(){
// Sets initial positions based on a crystal grid
  int ix=0, iy=0, iz=0, spacing=1; 

  while ((spacing*spacing*spacing)<amount) spacing++;

  for (int i=0; i<amount; i++){
    particlelist[i].mass = mass;
    particlelist[i].index = i;
    particlelist[i].box.x = 0;
    particlelist[i].box.y = 0;
    particlelist[i].box.z = 0;
    particlelist[i].position.x = ((ix+0.5)*region)/spacing;
    particlelist[i].position.y = ((iy+0.5)*region)/spacing;
    particlelist[i].position.z = ((iz+0.5)*region)/spacing;
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
void Particles::init_velocities(){
// Sets initial velocities and removes net momentum
  double momentum[3] = {0,0,0};
  double r;

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for private(r)
  #endif
    for (int i=0; i<amount; i++){
      r = rand_range(-0.5,0.5);
      particlelist[i].velocity.x = r;
      momentum[0] += r;
      r = rand_range(-0.5,0.5);
      particlelist[i].velocity.y = r;
      momentum[1] += r;
      r = rand_range(-0.5,0.5);
      particlelist[i].velocity.z = r;
      momentum[2] += r;
    }
  #pragma omp barrier
  for (int i=0; i<3; i++){
    momentum[i] = momentum[i]/amount;
  }
  for (int i=0; i<amount; i++){
    particlelist[i].velocity.x -= momentum[0];
    particlelist[i].velocity.y -= momentum[1];
    particlelist[i].velocity.z -= momentum[2];
  }
}
void Particles::update_force(int ** nlist){
// Calculates the force based on the neighborlist
  double r = region/2.;

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; i++){
      particlelist[i].force.x = 0;
      particlelist[i].force.y = 0;
      particlelist[i].force.z = 0;
    }
  #pragma omp barrier

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; i++){
      int neighbors = nlist[i][0];
      for (int j=1; j<=neighbors; j++){
        int index = nlist[i][j];
        Vec pos_i = particlelist[i].position;
        Vec pos_j = particlelist[index].position;
        double d = distance_periodic(pos_i, pos_j, region);
        double di = 1./d;
        double d3 = di*di*di;
        double dx = pos_i.x - pos_j.x;
        if (dx>r){
          dx -= region;
        }
        else if (dx<-r){
          dx += region;
        }
        double dy = pos_i.y - pos_j.y;
        if (dy>r){
          dy -= region;
        }
        else if (dy<-r){
          dy += region;
        }
        double dz = pos_i.z - pos_j.z;
        if (dz>r){
          dz -= region;
        }
        else if (dz<-r){
          dz += region;
        }
        particlelist[i].force.x += 48.*d3*(d3-0.5)*di*dx;
        particlelist[i].force.y += 48.*d3*(d3-0.5)*di*dy;
        particlelist[i].force.z += 48.*d3*(d3-0.5)*di*dz;
      }
    }
  #pragma omp barrier
}
void Particles::update_velocity(Vec ** xieta, double dt, double g, double sigma){
// Updates the velocity based on a second order integrator
// for the langevin equation
  double s3 = 1./sqrt(3);
  double sdt = sqrt(dt);
  double s3dt = sqrt(dt*dt*dt);

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; i++){
      double xi_x = xieta[i][0].x;
      double xi_y = xieta[i][0].y;
      double xi_z = xieta[i][0].z;
      double eta_x = xieta[i][1].x;
      double eta_y = xieta[i][1].y;
      double eta_z = xieta[i][1].z;
      double vel_x = particlelist[i].velocity.x;
      double vel_y = particlelist[i].velocity.y;
      double vel_z = particlelist[i].velocity.z;
      double force_x = particlelist[i].force.x;
      double force_y = particlelist[i].force.y;
      double force_z = particlelist[i].force.z;
      particlelist[i].velocity.x = vel_x + (0.5*dt*force_x)/mass - 0.5*dt*g*vel_x + 0.5*sdt*sigma*xi_x -
					0.125*dt*dt*g*(force_x/mass - g*vel_x) - 0.25*s3dt*g*sigma*(0.5*xi_x + s3*eta_x);
      particlelist[i].velocity.y = vel_y + (0.5*dt*force_y)/mass - 0.5*dt*g*vel_y + 0.5*sdt*sigma*xi_y -
					0.125*dt*dt*g*(force_y/mass - g*vel_y) - 0.25*s3dt*g*sigma*(0.5*xi_y + s3*eta_y);
      particlelist[i].velocity.z = vel_z + (0.5*dt*force_z)/mass - 0.5*dt*g*vel_z + 0.5*sdt*sigma*xi_z -
					0.125*dt*dt*g*(force_z/mass - g*vel_z) - 0.25*s3dt*g*sigma*(0.5*xi_z + s3*eta_z);
    }
  #pragma omp barrier
}
void Particles::update_position(Vec ** xieta, double dt, double sigma){
// Updates the velocity based on a second order integrator
// for the langevin equation
  double s3dt = sqrt(dt*dt*dt);
  double sq12 = 1./sqrt(12);

  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
    for (int i=0; i<amount; i++){
      double eta_x = xieta[i][1].x;
      double eta_y = xieta[i][1].y;
      double eta_z = xieta[i][1].z;
      double vel_x = particlelist[i].velocity.x;
      double vel_y = particlelist[i].velocity.y;
      double vel_z = particlelist[i].velocity.z;
      double new_x = particlelist[i].position.x;
      double new_y = particlelist[i].position.y;
      double new_z = particlelist[i].position.z;
      new_x += dt*vel_x + s3dt*sigma*sq12*eta_x;
      new_y += dt*vel_y + s3dt*sigma*sq12*eta_y;
      new_z += dt*vel_z + s3dt*sigma*sq12*eta_z;
      new_x = wrap_position(new_x);
      new_y = wrap_position(new_y);
      new_z = wrap_position(new_z);
      particlelist[i].position.x = new_x;
      particlelist[i].position.y = new_y;
      particlelist[i].position.z = new_z;
    }
  #pragma omp barrier
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
Particles::~Particles(void){
// Particle destructor - deallocates memory
  particlelist = 0;
  delete[] particlelist;
}

class Neighborlist{
  public:
    int ** nlist;
    void update_neighborlist(Atom * particlelist);
    void init(int amount_in, int refresh_in, int num_threads_in, double dt_in, double region_in, double radius_in);
    ~Neighborlist();
  private:
    int amount;
    int refresh;
    int called;
    int num_threads;
    int expected_neighbors;
    int * amount_neighbor;
    double radius;
    double region;
    double rv;
    double vmax;
    double dt;
    double max_abs(Atom * particlelist);
    int expected(double radius);
    int max_neighbors(void);
    void update_radius(Atom * particlelist);
    void allocate_mem(int exp_neighbors);
};
void Neighborlist::init(int amount_in, int refresh_in, int num_threads_in, double dt_in, double region_in, double radius_in){
  amount = amount_in;
  refresh = refresh_in;
  num_threads = num_threads_in;
  dt = dt_in;
  region = region_in;
  radius = radius_in;
  called = 0;
  amount_neighbor = new int[amount];
  expected_neighbors = expected(radius);
  //allocate_mem(expected_neighbors);
  allocate_mem(amount);
}
int Neighborlist::expected(double radius_in){
// Calculates the expected amount of neighbors and adds 30%
  double r  = radius_in;
  int exp_neigh = (((4./3.)*PI*r*r*r)/(region*region*region))*amount;
  exp_neigh = exp_neigh*1.3;
  return exp_neigh;
}
void Neighborlist::allocate_mem(int exp_neighbors){
// Allocates memory for the neighborlist
  nlist = new int*[amount];
  for (int i=0; i<amount; ++i){
    nlist[i] = new int[exp_neighbors];
  }
}
void Neighborlist::update_neighborlist(Atom * particlelist){
// Updates the neighborlist
  if (called%refresh==0){
    update_radius(particlelist);
    int max = max_neighbors();
    if (max >= 0.8*expected_neighbors){
      //allocate_mem(max);
      expected_neighbors = max;
    }
  #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #pragma omp parallel for
  #endif
      for (int i=0; i<amount; i++){
        int l = 0;
        for (int j=0; j<amount; j++){
          if (i!=j){
            double distance_ij = distance_periodic(particlelist[i].position, particlelist[j].position, region);
            if (distance_ij<rv*rv){
              nlist[i][l+1] = particlelist[j].index;
              l++;
            }
          }
        }
        amount_neighbor[i] = l;
        nlist[i][0] = amount_neighbor[i];
      }
    #pragma omp barrier
  }  
  called++;
}
void Neighborlist::update_radius(Atom * particlelist){
// Updates the cutoff radius based on the maximum absolute velocity
  vmax = max_abs(particlelist);
  rv = radius + refresh*vmax*dt;
}
int Neighborlist::max_neighbors(){
// Returns the maximum amount of neighbors for the particles
  int max1, max2=0;
  for (int i=0; i<amount; i++){
    max1 = amount_neighbor[i];
    if (max1 > max2){
      max2 = max1;
    }
  }
  return max2;
}
double Neighborlist::max_abs(Atom * particlelist){
// Calculates the maximum abolute value squared of all velocities 
  double c, max1=0, max2=0;
  
  for (int i=0; i<amount; i++){
    c = particlelist[i].velocity.x;
    max1 += c*c;
    c = particlelist[i].velocity.y;
    max1 += c*c;
    c = particlelist[i].velocity.y;
    max1 += c*c;
    if (max1 >= max2){
      max2 = max1;
      max1 = 0;
    }
    else {
      max1 = 0;
    }
  }
  return max2;
}
Neighborlist::~Neighborlist(void){
// Neighborlist destructor - deallocates memory
  for (int i=0; i<amount; ++i){
    nlist[i] = 0;
    delete[] nlist[i];
  }
  nlist = 0;
  delete[] nlist;
}

class System{
  public:
    int amount;
    int refresh;
    int timesteps;
    int num_threads;
    double dt;
    double kB;
    double T;
    double m;
    double rho;
    double radius;
    double region;
    double g;
    double sigma;
    double seed;
    Particles particles;
    Neighborlist neighbor;
    XiEta rnd_list;
    System(int amount_in, int refresh_in, int timesteps_in, int num_threads_in, double dt_in, double kB_in, double T_in, double m_in, double rho_in, double radius_in, double g_in, double seed_in);
    void update_neighborlist();
    void update_force();
    void update_velocity();
    void update_position();
    void calc_rnd();
};
System::System(int amount_in, int refresh_in, int timesteps_in, int num_threads_in, double dt_in, double kB_in, double T_in, double m_in, double rho_in, double radius_in, double g_in, double seed_in){
  amount = amount_in;
  refresh = refresh_in;
  timesteps = timesteps_in;
  num_threads = num_threads_in;
  dt = dt_in;
  kB = kB_in;
  T = T_in;
  m = m_in;
  rho = rho_in;
  radius = radius_in;
  g = g_in;
  seed = seed_in;
  srand(seed);
  region = pow((amount*m)/rho,1./3.);
  sigma = sqrt((2*kB*T*g)/m);
  particles.init(amount,num_threads,region,m);
  particles.init_positions();
  particles.init_velocities();
  neighbor.init(amount,refresh,num_threads,dt,region,radius);
  rnd_list.init(amount,num_threads);
}
void System::update_neighborlist(void){
  neighbor.update_neighborlist(particles.particlelist);
}
void System::update_force(void){
  particles.update_force(neighbor.nlist);
}
void System::calc_rnd(void){
  rnd_list.calc();
}
void System::update_velocity(void){
  particles.update_velocity(rnd_list.xieta,dt,g,sigma);
}
void System::update_position(void){
  particles.update_position(rnd_list.xieta,dt,sigma);
}
