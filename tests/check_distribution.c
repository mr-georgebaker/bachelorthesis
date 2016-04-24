// Check if function returns normally distributed values
// and creates a file "data.dat" which can be used for plotting
// with gnuplot

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.14159265

/////////////
//Datatypes//
/////////////

// Datatype for storing two normaly distributed values
typedef struct
{
  double xi; 
  double eta;
} XiEta;

////////////////////
//Helper functions//
////////////////////

double rand_range(double min, double max){
// Returns a random number between min and max
  double scaled = (double)rand()/RAND_MAX;
  return scaled*(max-min)+min;
}

XiEta normal_value(){
// Returns two normaly distributed, independent values with mu=0 and sigma=1
// based on the Box-Muller Method
  double u = rand_range(0,1);
  double v = rand_range(0,1);

  XiEta x = {sqrt(-2*log(u))*cos(2*PI*v), sqrt(-2*log(u))*sin(2*PI*v)};
  return x;
}

int main(){
  //////////////////
  //Initialization//
  //////////////////
  
  int i;
  int dim = 100000;
  XiEta rnd_dis[dim];
  FILE *fp;
  fp = fopen("data.dat", "w+");

  ///////////////
  //Calculation//
  ///////////////

  for (i=0; i<dim; i++){
    rnd_dis[i] = normal_value();
    fprintf(fp, "%f\n",rnd_dis[i].xi);
  }
  
  fclose(fp);

  return 0;
}
