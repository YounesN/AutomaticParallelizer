#include <iostream>
#include "stencil.hpp"
#include <ctime>

int main(int argc, char * argv[])
{
  // Check the number of parameters
  if(argc < 3) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] <<
      " <filename> <iterations>" << std::endl;
    return -1;
  }

  // store the filename given by user into a string
  std::string filename(argv[1]);

  // get the number of iterations or time
  int iterations = atoi(argv[2]);

  // Create an object of stencil
  Stencil stencil(filename, iterations);

  // read data from the file
  stencil.ReadData();

  // run stencil code for 20 iterations
  std::clock_t c_start = std::clock();
  stencil.RunStencil();
  std::clock_t c_end = std::clock();

  // print time took to run the algorithm
  long double time_elapsed_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
  std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";

  // output data to another file
  std::string output_name(argv[1]);
  output_name = "output_" + output_name;
  stencil.OutputData(output_name);

  return 0;
}
