/*
	Molecular dynamics simulation using a second order integrator and neighbor list
	Creates a file positions.xyz after the simulation containing the positions of
	the particles to use a molecular editor to plot them in the form

	Amount of particles
	Description
	Element  x-position	y-position	z-position
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 8
#define PI 3.14159265

/////////////
//Datatypes//
/////////////

typedef struct
// Datatype for storing two normaly distributed values
{
  double xi; 
  double eta;
} XiEta;

////////////////////
//Helper functions//
////////////////////

double rand_range(double min, double max){
// Returns a uniformally distributed random number between min and max
  double scaled = (double)rand()/RAND_MAX;
  return scaled*(max-min)+min;
}

XiEta normal_value(double mu, double sigma){
// Returns two normaly distributed, independent values with mu=0 and sigma=1
// based on the Box-Muller Method
  double u = rand_range(0,1);
  double v = rand_range(0,1);

  XiEta x = {mu + sigma*sqrt(-2*log(u))*cos(2*PI*v), mu + sigma*sqrt(-2*log(u))*sin(2*PI*v)};
  return x;
}

double distance_periodic(double *xi, double *xj, double region){
// Returns the squared distance between two points in three dimensions
// for periodic boundary conditions
  double d=0, r=region/2., c;
  int i;

  for (i=0; i<3; i++){
    c = xi[i]-xj[i];
    if (c > r){
      c -= region;
    }
    else if (c < -r){
      c += region;
    }
    d += c*c;
  }
  return d;
}

double *sum_vel(double velocity[][3], int amount){
// Returns the sum of all velocities (should be 0 if g=0)
  int i,j;
  static double vel_sum[3] = {0,0,0};

  for (i=0; i<3; i++){
    for (j=0; j<amount; j++){
      vel_sum[i] += velocity[j][i];
    }
  }
  return vel_sum;
}

double max_abs(double velocities[][3], int amount){
// Returns the maximum absolute squared of an array
  int i, j;
  double max1 = 0;
  double max2 = 0;
  double c;
 
  for (i=0; i<amount; i++){
    for (j=0; j<3; j++){
      c = velocities[i][j];
      max1 += c*c;
    }
    //max1 = sqrt(max1);
    if (max1 >= max2){
      max2 = max1;
      max1 = 0;
    }
    else{
      max1 = 0;
    }
  }
  return max2; 
}

/////////////////////////
//Measurement functions//
/////////////////////////

double kinetic_energy(double velocities[][3], double m, int amount){
// Returns the kinetic energy
  int i,j;
  double vel, kinetic = 0;

  for (i=0; i<amount; i++){
    for (j=0; j<3; j++){
      vel = velocities[i][j];
      kinetic += (m*vel*vel)/2.;
    }
  }

  return kinetic;
}

///////////////////
//Setup functions//
///////////////////

void create_particles(double positions[][3], int amount, double region){
// Creates particles equally spaced inside the box
  int iix=0, iiy=0, iiz=0, n3=1, i;

  while ((n3*n3*n3)<amount) n3++;

  for (i=0;i<amount;i++) {
    positions[i][0] = ((double)iix+0.5)*region/n3;
    positions[i][1] = ((double)iiy+0.5)*region/n3;
    positions[i][2] = ((double)iiz+0.5)*region/n3;
    iix++;
    if (iix==n3) {
      iix=0;
      iiy++;
      if (iiy==n3) {
	iiy=0;
	iiz++;
      }
    }
  }
}

void initial_velocities(double velocities[][3], double sigma, double m, int amount){
// Sets the initial velocity based on the temperature and removes net momentum
  int i, j;
  double momentum[3] = {0,0,0};

  for (i=0; i<amount; i++){
    for (j=0; j<3; j++){
      velocities[i][j] = normal_value(0,sigma).xi;
      momentum[j] += velocities[i][j]*m;
    }
  }
  for (j=0; j<3; j++){
    momentum[j] = momentum[j]/((double)amount*m);
  }
  for (i=0; i<amount; i++){
    for (j=0; j<3; j++){
      velocities[i][j] -= momentum[j];
    }
  }
}

////////////////////
//Update functions//
////////////////////

void update_neighborlist(double ****neighbor, double *amount_neighbor, double particles[][3], double radius, int amount, double region){
// Updates the neighbor list for each particle within a given radius
// The first element of each row contains a pointer to the amount of neighbors
  int i,j,k;

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j,k)
    for (i=0; i<amount; i++){
      int l = 0;
      for (j=0; j<amount; j++){
        if (i!=j){
          double distance_ij = distance_periodic(particles[i],particles[j],region);
          if (distance_ij<radius*radius){
            for (k=0; k<3; k++){
              neighbor[i][l+1][k] = &particles[j][k];
            }
            l++;
          }
        }
      }
      amount_neighbor[i] = (double)l;
      neighbor[i][0][0] = &amount_neighbor[i];
    }
  #pragma omp barrier
}

void calculate_force(double forces[][3], double ****neighbor, double particles[][3], double region, int amount){
// Updates the force for each particle based on the neighbor list
  int i,j,k;
  double r = region/2.;

  
  for (i=0; i<amount; i++){
    for (j=0; j<3; j++){
      forces[i][j] = 0;
    }
  }
  
  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j,k)
    for (i=0; i<amount; i++){
      int neighbors = *neighbor[i][0][0];
      for (j=1; j<=neighbors; j++){
        double d = distance_periodic(particles[i],*neighbor[i][j],region);
        double di = 1./d;
        double d3 = di * di *di;
        for (k=0; k<3; k++){
          double dl = particles[i][k] - *neighbor[i][j][k];
          if (dl > r){
            dl -= region;
          }
          else if (dl < -r){
            dl += region;
          }
          forces[i][k] += 48.*d3*(d3-0.5)*di*dl;
        }
      }
    }
  #pragma omp barrier
}

void update_velocity(double velocities[][3], double forces[][3], double dt, double g, double sigma, double m, XiEta rnd_XiEta[][3], int amount){
// Updates the velocities
  int i,j;
  
  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j)
    for (i=0; i<amount; i++){
      for (j=0; j<3; j++){
        double vel = velocities[i][j];
        double xi = rnd_XiEta[i][j].xi;
        double eta = rnd_XiEta[i][j].eta;
        double force = forces[i][j];
        velocities[i][j] = vel + 0.5*dt*force/m - 0.5*dt*g*vel + 0.5*sqrt(dt)*sigma*xi -
				0.125*dt*dt*g*(force/m - g*vel) - 0.25*ldexp(dt,3./2.)*g*sigma*(0.5*xi + (1./sqrt(3.))*eta);
      }
    }
  #pragma omp barrier  
}

void update_position(double positions[][3], double velocities[][3], double dt, double sigma, double region, XiEta rnd_XiEta[][3], int amount){
// Updates the positions
  int i,j;

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j)
    for (i=0; i<amount; i++){
      for (j=0; j<3; j++){
        double actual = positions[i][j];
        actual += dt*velocities[i][j] + ldexp(dt,3./2.)*sigma*(1./sqrt(12.))*rnd_XiEta[i][j].eta;
        if (actual < -0.5*region){
          positions[i][j] = actual + region;
        }
        else if (actual >= 0.5*region){
          positions[i][j] = actual - region;
        }
        else {
          positions[i][j] = actual;
        }
      }
    }
  #pragma omp barrier  
}

////////
//MAIN//
////////

int main(){

  ///////////////////
  //Initialize data//
  ///////////////////

  time_t t;
  int i,j,k;
  int amount = 1000;				// Amount of particles
  int timesteps = 1000;				// Amount of timesteps
  int refresh = 10;				// Amount of timesteps without updating neighbour list
  double rc = 2.5; 				// Cutoff radius
  double rv;					// Neighborlist radius
  double rho = 0.5;				// Density
  double m = 1;				  	// Mass of particles
  double region = pow((amount*m)/rho,1./3.);	// Sidelength of region
  double dt = 0.005;				// Timestep
  double kB = 1;				// Boltzmannconstant
  double g = 0;					// Guiding factor
  double T = 1;					// Temperature
  double sigma = sqrt((2*kB*T*g)/m);		// Whitenoise factor
  double sig = sqrt((kB*T)/m);			// Initial velocity factor
  double vmax;					// Maximum velocity of particles
  double kinetic;				// Kinetic energy
  double positions[amount][3];			// Particle positions
  double velocities[amount][3];			// Particle velocities
  double forces[amount][3];			// Particle forces
  double amount_neighbors[amount];		// Amount of neighbors per particle
  double ****neighbor;				// Neighbor List
  double *vel_sum;				// For checking momentum conservation
  XiEta rnd_XiEta[amount][3];			// Array of linear independent normally distributed numbers
  FILE *fp;

  srand((unsigned) time(&t));			// Set seed for random number generator
  fp = fopen("positions.xyz", "w+");

  ///////////////////////////
  //Setup initial positions//
  ///////////////////////////
  
  create_particles(positions, amount, region);
  initial_velocities(velocities, sig, m, amount);
  
  ///////////////////////////////
  //Allocate memory on the heap//
  ///////////////////////////////

  neighbor = malloc(amount*sizeof(double ***));
  if (neighbor){
    for (i=0; i<amount; ++i){
      neighbor[i] = malloc(amount*sizeof(double **));
      if (neighbor[i]){
        for (j=0; j<amount; ++j){
          neighbor[i][j] = malloc(3*sizeof(double *));
          if (!neighbor[i][j]){
            printf("\nMemory allocation error!\n");
          }
        }
      }
    }
  }

  //////////////////
  //Main algorithm//
  //////////////////

  for (i=0; i<timesteps; i++){
    if (i%refresh==0){
      vmax = sqrt(max_abs(velocities, amount));
      rv = rc + refresh*vmax*dt;
      update_neighborlist(neighbor, amount_neighbors, positions, rv, amount, region);  
    }
    for (j=0; j<amount; j++){
      for (k=0; k<3; k++){
        rnd_XiEta[j][k] = normal_value(0,1);
      } 
    }

    calculate_force(forces, neighbor, positions, region, amount);
    update_velocity(velocities, forces, dt, g, sigma, m, rnd_XiEta, amount);
    update_position(positions, velocities, dt, sigma, rnd_XiEta, amount);
    calculate_force(forces, neighbor, positions, region, amount);
    update_velocity(velocities, forces, dt, g, sigma, m, rnd_XiEta, amount);
  }

fprintf(fp, "%d\n", amount);
fprintf(fp, "Argon fluid positions\n");
  
for (i=0; i<amount; i++){
  fprintf(fp, "Ar\t");
  for (j=0; j<3; j++){
    fprintf(fp, "%f\t", positions[i][j]);
  }
  fprintf(fp, "\n"); 
}

close(fp);

  

  ////////////////////////////
  //Clear memory on the heap//
  ////////////////////////////

  for (i=0; i<amount; i++){
    for (j=0; j<amount; j++){
      free(neighbor[i][j]);
    }
    free(neighbor[i]);
  }
  free(neighbor);

  return 0;
}



