#include <iostream>
#include "stencil.hpp"
#include <ctime>
#include <omp.h>
#include "gettime.h"

int main(int argc, char * argv[])
{
  // Check the number of parameters
  if(argc < 4) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <filename> <iterations> <number_of_threads>" << std::endl;
    return -1;
  }

  // store the filename given by user into a string
  std::string filename(argv[1]);

  // get the number of iterations or time
  int iterations = atoi(argv[2]);

#ifdef _OPENMP
    // set 4 threads for now
  int numThreads = atoi(argv[3]);
  omp_set_num_threads(numThreads);

  // print some OpenMP info
  std::cout << "OpenMP is activated!" << std::endl;
  std::cout << "Number of threads: " << numThreads << std::endl;
#else
    std::cout << "OpenMP is not activated!" << std::endl;
#endif

  // Create an object of stencil
  Stencil stencil(filename, iterations);

  // read data from the file
  stencil.ReadData();

  // run stencil code for 20 iterations
  double wall0 = get_wall_time();
  double cpu0  = get_cpu_time();
  stencil.RunStencil();
  double wall1 = get_wall_time();
  double cpu1  = get_cpu_time();

  // print time took to run the algorithm
  std::cout << "Wall time: " << wall1 - wall0 << " s\n";
  std::cout << "CPU time: " << cpu1 - cpu0 << " s\n";

  // output data to another file
  std::string output_name(argv[1]);
  output_name = "output_" + output_name;
  stencil.OutputData(output_name);

  return 0;
}