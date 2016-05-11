#include "datastructure.h"
#include <fstream>
#include <getopt.h>
using namespace std;

void print_usage() {
    std::cout << "Usage:\t-a \t Amount of particles \n\t\
-t \t Amount of timesteps\n\t\
-r \t Neighborlist refresh rate\n\t\
-d \t Stepsize\n\t\
-T \t Temperature\n\t\
-m \t Mass\n\t\
-D \t Density\n\t\
-c \t Cutoff-radius\n\t\
-g \t Whitenoise factor\n\t\
-s \t Seed\n\t\
-o \t Number of threads (optional)\n\t\
-f \t Filename for positions (optional)" << endl;
}

int main(int argc, char *argv[]){

  double kB = 1;
  int opt = 0;
  int amount = -1;
  int refresh;
  int timesteps;
  int num_threads;
  double dt;
  double T;
  double m;
  double rho;
  double radius;
  double g;
  double seed;
  char* filename_pos;
  bool pos_file = false;

  static struct option long_options[] = {
        {"amount",   	required_argument, 0,  'a' },
        {"timesteps", 	required_argument, 0,  't' },
        {"refresh",   	required_argument, 0,  'r' },
        {"dt",   	required_argument, 0,  'd' },
        {"T",   	required_argument, 0,  'T' },
        {"m",   	required_argument, 0,  'm' },
        {"rho",   	required_argument, 0,  'D' },
        {"r",   	required_argument, 0,  'c' },
        {"g",   	required_argument, 0,  'g' },
        {"seed",   	required_argument, 0,  's' },
        {"threads", 	optional_argument, 0,  'o' },
        {"file",   	optional_argument, 0,  'f' },   
    };

  int long_index =0;
  while ((opt = getopt_long(argc, argv,"a:t:r:d:T:m:D:c:g:s:o:f:M:", 
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
      case 'T' : T = atof(optarg);
        break;
      case 'm' : m = atof(optarg);
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
      default: print_usage(); 
        exit(EXIT_FAILURE);
      }
    }
    if (amount < 0) {
      print_usage();
      exit(EXIT_FAILURE);
    }
  

  ofstream file_pos;
  file_pos.open(filename_pos);
  
  System system(amount,refresh,timesteps,num_threads,dt,kB,T,m,rho,radius,g,seed);

  double kin, avk=0, avk2=0, pot;
  Vec momentum;

  for (int i=1; i<=timesteps; ++i){
    if (pos_file){
      file_pos << amount << endl;
      file_pos << "Argon fluid simulation" << endl;
      for (int j=0; j<amount; j++){
        file_pos << "Ar" << "\t" << system.particles.particlelist[j].position.x << "\t" << system.particles.particlelist[j].position.y << "\t" << system.particles.particlelist[j].position.z << "\n"; 
      }
    }
    system.update_neighborlist();
    system.update_force();
    system.calc_rnd();
    system.update_velocity();
    system.update_position();
    system.update_force();
    system.update_velocity();
  }

  
  file_pos.close();
  return 0;
}
