/*
Creates a neighbor list for each particle
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
  int i,j,k,l;
  int p = 500; // Amount of particles
  double random[p][dim], *neighbor[p][p-1][dim];
  double r = 20; // Cutoff radius
  double distance_ij;

  // Create random particles
  for (i=0; i<p; i++){
    for (j=0; j<dim; j++){
      random[i][j] = rand() % 100;
    }
  }

  // Create neighbor list with pointers to positions
  for (i=0; i<p; i++){
    l = 0;
    for (j=0; j<p; j++){
      if (i!=j){
        distance_ij = distance(random[i],random[j],dim);
        if (distance_ij<r){
          for (k=0; k<dim; k++){
            neighbor[i][l][k] = &random[j][k];
          }
          l++;
        }
      }
    }
  }
 
  return 0;
}
