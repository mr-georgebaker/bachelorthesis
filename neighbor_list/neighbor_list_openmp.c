/*
Creates a neighbor list for all particles using parallelization with OpenMP
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 4

double distance(double xi[], double xj[], int l){
// Returns the (euclidian) distance between two points
  double d=0;
  int i;
  for (i = 0; i<l; i++){
    d += pow(xi[i]-xj[i],2);
  }
  d = sqrt(d);
  return d;
}

int main() 
{
  int dim = 3; // Dimension
  int i,j,k;
  int p = 500; // Amount of particles
  double particles[p][dim], *neighbor[p][p-1][dim];
  double r = 20; // Cutoff radius
  double distance_ij;

  // Create random particles
  for (i=0; i<p; i++){
    for (j=0; j<dim; j++){
      particles[i][j] = rand() % 100;
    }
  }

  // Create neighbor list with pointers to positions (OpenMP)
  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j,k)
    for (i=0; i<p; i++){
      int l = 0;
      for (j=0; j<p; j++){
        if (i!=j){
          distance_ij = distance(particles[i],particles[j],dim);
          if (distance_ij<r){
            for (k=0; k<dim; k++){
              neighbor[i][l][k] = &particles[j][k];
            }
            l++;
          }
        }
      }
    }
 
  return 0;
}