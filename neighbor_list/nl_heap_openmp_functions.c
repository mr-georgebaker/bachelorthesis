/*
Creates a neighbor list on the heap for all particles using parallelization with OpenMP
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#define NUM_THREADS 4

double distance_periodic(double * xi, double * xj, double region){
// Returns the distance between two points in three dimensions
// for periodic boundary conditions
  double d=0, c;
  int i;
  for (i=0; i<3; i++){
    c = xi[i]-xj[i];
    d += c*c;
  }
  d = sqrt(d);
  if (d > region/2){
    d -= region/2;
  }
  return d;
}

double create_neighbor_list(double ****neighbor, double *amount_neighbor, double particles[][3], double radius, int amount, double region){
// Returns the neighbor list for each particle within a given radius
// The first element of each row contains a pointer to the amount of neighbors
  int i,j,k;
  double distance_ij;

  omp_set_num_threads(NUM_THREADS);
  #pragma omp parallel for private(j,k)
    for (i=0; i<amount; i++){
      int l = 0;
      for (j=0; j<amount; j++){
        if (i!=j){
          distance_ij = distance_periodic(particles[i],particles[j],region);
          if (distance_ij<radius){
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

int main() 
{
  int i,j;
  int p = 1000; // Amount of particles
  double particles[p][3];
  double r = 20; // Cutoff radius
  double region = 100;
  double ****neighbor;
  double *amount_neighbor;

  amount_neighbor = malloc(p*sizeof(double));

  neighbor = malloc(p*sizeof(double ***));
  if (neighbor){
    for (i=0; i<p; ++i){
      neighbor[i] = malloc(p*sizeof(double **));
      if (neighbor[i]){
        for (j=0; j<p; ++j){
          neighbor[i][j] = malloc(3*sizeof(double *));
          if (!neighbor[i][j]){
            printf("\nMemory allocation error!\n");
          }
        }
      }
    }
  }

  // Create random particles
  for (i=0; i<p; i++){
    for (j=0; j<3; j++){
      particles[i][j] = rand() % (int)region;
    }
  }

  create_neighbor_list(neighbor, amount_neighbor, particles, r, p, region);
  
  // Free memory of array on the heap
  for (i=0; i<p; i++){
    for (j=0; j<p; j++){
      free(neighbor[i][j]);
    }
    free(neighbor[i]);
  }
  free(neighbor);
  
  free(amount_neighbor);

  return 0;
}
