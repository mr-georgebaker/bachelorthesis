#include "datastructure2.h"
#include <getopt.h>

void print_usage() {
    std::cerr << "Usage:\t-a \t Amount of particles \n\t\
-t \t Amount of timesteps (Default: 1)\n\t\
-r \t Neighborlist refresh rate (Default: 1)\n\t\
-d \t Stepsize (Default: 0.001)\n\t\
-T \t Temperature (Default: 1)\n\t\
-m \t Mass (Default: 1)\n\t\
-D \t Density (Default: 0.6)\n\t\
-c \t Cutoff-radius (Default: boxsize/2)\n\t\
-g \t g factor (Default: 0)\n\t\
-s \t Seed (Default: 1)\n\t\
-o \t Number of threads (optional) (Default: 1)\n\t\
-f \t Filename for position output (optional)\n\t\
-p \t Filename for position input (optional)\n\t\
-v \t Filename for velocity input (optional)" << std::endl;
}

int main(int argc, char *argv[]){

  int opt          = 0;
  int amount       = -1;
  int refresh      = 1;
  int timesteps    = 1;
  int num_threads  = 1;
  double dt        = 0.001;
  double temp      = 1.0;
  double mass      = 1.0;
  double rho       = 0.6;
  double radius    = 0;
  double g  	   = 0; 
  double seed      = 1;
  char* filename_pos;
  bool pos_file = false;
  bool init_pos = false;
  bool init_vel = false;
  char* init_pos_file;
  char* init_vel_file;

  static struct option long_options[] = {
        {"amount",   	required_argument, NULL,  'a' },
        {"timesteps", 	required_argument, NULL,  't' },
        {"refresh",   	required_argument, NULL,  'r' },
        {"dt",   	required_argument, NULL,  'd' },
        {"T",   	required_argument, NULL,  'T' },
        {"m",   	required_argument, NULL,  'm' },
        {"rho",   	required_argument, NULL,  'D' },
        {"r",   	required_argument, NULL,  'c' },
        {"g",   	required_argument, NULL,  'g' },
        {"seed",   	required_argument, NULL,  's' },
        {"threads", 	optional_argument, NULL,  'o' },
        {"file",   	optional_argument, NULL,  'f' },   
        {"pos_file",   	optional_argument, NULL,  'p' }, 
        {"vel_file",   	optional_argument, NULL,  'v' }, 
    };

  int long_index =0;
  while ((opt = getopt_long(argc, argv,"a:t:r:d:T:m:D:c:g:s:o:f:M:p:v:", 
                   long_options, &long_index )) != -1) {
    switch (opt) {
      case 'a' : amount = atoi(optarg);
        break;
      case 't' : timesteps = atoi(optarg);
        break;
      case 'r' : refresh = atoi(optarg);
        break;
      case 'd' : dt = atof(optarg);
        break;
      case 'T' : temp = atof(optarg);
        break;
      case 'm' : mass = atof(optarg);
        break;
      case 'D' : rho = atof(optarg);
        break;
      case 'c' : radius = atof(optarg);
        break;
      case 'g' : g = atof(optarg);
        break;
      case 's' : seed = atof(optarg);
        break;
      case 'f' : filename_pos = optarg;
        pos_file = true;
        break;
      case 'o' : num_threads = atoi(optarg);
        break;
      case 'p' : init_pos_file = optarg;
        init_pos = true;
        break;
      case 'v' : init_vel_file = optarg;
        init_vel = true;
        break;
      default: print_usage(); 
        exit(EXIT_FAILURE);
      }
    }
    if (amount < 0) {
      print_usage();
      exit(EXIT_FAILURE);
    }

  //////////
  // MAIN //
  //////////

  System system(amount,refresh,timesteps,num_threads,dt,temp,mass,rho,radius,g,seed,init_pos,init_pos_file,init_vel,init_vel_file);
  system.print_init();
  system.update_neighborlist();
  system.update_forces();

  // Initial energy
  double E0 = system.potential_energy() + system.kinetic_energy();

  for (int i=0; i<timesteps; ++i){
    //system.print_positions();
    //system.print_forces();
    //system.print_velocities();
    //system.print_energy();
    //system.print_total_momentum();
    //system.print_neighborlist();

    if (g!=0){
      system.update_rnd();
    }
    system.update_velocities();
    system.update_positions();
    system.update_forces();
    system.update_velocities();
    system.update_neighborlist();

    std::cout << i << "\t" << system.pressure() << std::endl;
  }

  return 0;
}
