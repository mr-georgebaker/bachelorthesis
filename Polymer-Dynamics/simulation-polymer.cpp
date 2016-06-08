/*
	Simulation of a polymer using the bead-spring model,
	the Lennard-Jones-Potential and an external potential.
	Integrator by Eric Vanden-Eijnden and Giovanni Ciccotti
	for the Langevin-Equation.
*/

#include "datastructure-polymer.h"
#include <getopt.h>

void print_usage() {
    std::cerr << "Usage:\t-a \t Amount of particles\n\t\
-t \t Amount of timesteps (Default: 1)\n\t\
-h \t Stepsize (Default: 0.001)\n\t\
-T \t Temperature (Default: 1)\n\t\
-m \t Mass (Default: 1)\n\t\
-g \t Friction factor (Default: 0)\n\t\
-d \t Average bond length (Default: 2)\n\t\
-b \t Boxsize (Default: 20)\n\t\
-r \t Neighborlist refresh rate (Default: 1)\n\t\
-c \t Cutoff-radius (Default: boxsize/2)\n\t\
-s \t Seed (Default: 1)\n\t\
-e \t External potential width (Default: 0)\n\t\
-z \t External potential root (Default: 2)\n\t\
-o \t Number of threads (optional) (Default: 1)\n\t\
-f \t Filename for position output (optional)\n\t\
-p \t Filename for position input (optional)\n\t\
-v \t Filename for velocity input (optional)" << std::endl;
}

int main(int argc, char *argv[]){

  int opt		= 0;
  int amount		= -1;
  int timesteps		= 1;
  int refresh		= 1;
  int threads		= 1;

  double dt 		= 0.001;
  double temp		= 1.;
  double mass		= 1.;
  double friction	= 0.;
  double distance	= 2.;
  double region_length	= 20.0;
  double radius_c	= region_length/2.;
  double seed		= 1.;
  double epsilon	= 0.;
  double c		= 2.;

  char * filename_pos_out;
  char * filename_pos_in;
  char * filename_vel_in;

  bool file_pos_out	= false;
  bool file_pos_in	= false;
  bool file_vel_in	= false;

  static struct option long_options[] = {
    {"amount",		required_argument,	NULL,	'a'},
    {"timesteps",	required_argument,	NULL,	't'},
    {"stepsize",	required_argument,	NULL,	'h'},
    {"temperature",	required_argument,	NULL,	'T'},
    {"mass",		required_argument,	NULL,	'm'},
    {"friction",	required_argument,	NULL,	'g'},
    {"distance",	required_argument,	NULL,	'd'},
    {"boxsize",		required_argument,	NULL,	'b'},
    {"refresh-rate",	required_argument,	NULL,	'r'},
    {"cutoff-radius",	required_argument,	NULL,	'm'},
    {"seed",		required_argument,	NULL,	's'},
    {"epsilon",		required_argument,	NULL,	'e'},
    {"c",		required_argument,	NULL,	'z'},
    {"threads",		optional_argument,	NULL,	'o'},
    {"pos_out",		optional_argument,	NULL,	'f'},
    {"pos_in",		optional_argument,	NULL,	'p'},
    {"vel_in",		optional_argument,	NULL,	'v'},
  };

  int long_index =0;
  while ((opt = getopt_long(argc, argv,"a:t:h:T:m:g:d:b:r:m:s:e:z:o:f:p:v:", 
                   long_options, &long_index )) != -1) {
    switch (opt) {
      case 'a' : amount = atoi(optarg);
        break;
      case 't' : timesteps = atoi(optarg);
        break;
      case 'h' : dt = atof(optarg);
        break;
      case 'T' : temp = atof(optarg);
        break;
      case 'm' : mass = atof(optarg);
        break;
      case 'g' : friction = atof(optarg);
        break;
      case 'd' : distance = atof(optarg);
        break;
      case 'b' : region_length = atof(optarg);
        break;
      case 'r' : refresh = atoi(optarg);
        break;
      case 'c' : radius_c = atof(optarg);
        break;
      case 's' : seed = atof(optarg);
        break;
      case 'e' : epsilon = atof(optarg);
        break;
      case 'z' : c = atof(optarg);
        break;
      case 'o' : threads = atoi(optarg);
        break;
      case 'f' : filename_pos_out = optarg;
        file_pos_out = true;
        break;
      case 'p' : filename_pos_in = optarg;
        file_pos_in = true;
        break;
      case 'v' : filename_vel_in = optarg;
        file_vel_in = true;
        break;
      default: print_usage(); 
        exit(EXIT_FAILURE);
    }
  }
  if (amount < 0) {
    print_usage();
    exit(EXIT_FAILURE);
  }

  System system(amount,timesteps,refresh,threads,dt,temp,mass,friction,distance,region_length,radius_c,seed,epsilon,c,file_pos_out,filename_pos_out, file_pos_in,filename_pos_in,file_vel_in,filename_vel_in);
  system.print_init();

  system.print_positions();
  //system.print_velocities();

  system.update_neighborlist();
  system.update_forces();

  //system.print_end_to_end();
  //system.print_forces();
  //system.print_neighborlist();

  for (int i=0; i<timesteps; ++i){
    if (friction!=0){
      system.update_rnd();
    }
    system.update_velocities();
    system.update_positions();
    system.update_forces();
    system.update_velocities();
    system.update_neighborlist();

    system.print_positions();
    //system.print_forces();
    //system.print_velocities();
    //system.print_end_to_end();
  }
}
