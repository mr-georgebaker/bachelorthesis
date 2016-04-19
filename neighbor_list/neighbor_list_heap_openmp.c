/*
Creates a neighbor list on the heap for all particles using parallelization with OpenMP
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
  int p = 1000; // Amount of particles
  double particles[p][dim];
  double r = 20; // Cutoff radius
  double distance_ij;

  // Create array for neighbor list on the heap
  double ****neighbor;
  neighbor = malloc(p*sizeof(double ***));
  if (neighbor){
    for (i=0; i<p; ++i){
      neighbor[i] = malloc((p-1)*sizeof(double **));
      if (neighbor[i]){
        for (j=0; j<p-1; ++j){
          neighbor[i][j] = malloc(dim*sizeof(double *));
          if (!neighbor[i][j]){
            printf("\nMemory allocation error!\n");
          }
        }
      }
    }
  }

  // Create random particles
  for (i=0; i<p; i++){
    for (j=0; j<dim; j++){
      particles[i][j] = rand() % 100;
    }
  }

  // Create neighbor list with pointers to positions (OpenMP implementation)
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

  // Free memory of array on the heap
  for (i=0; i<p; i++){
    for (j=0; j<p-1; j++){
      free(neighbor[i][j]);
    }
    free(neighbor[i]);
  }
  free(neighbor);
 
  return 0;
}
