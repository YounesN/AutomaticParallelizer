#include <iostream>
#include "stencil.hpp"

int main(int argc, char * argv[])
{
  // Check the number of parameters
  if(argc < 2) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
    return -1;
  }

  // store the filename given by user into a string
  std::string filename(argv[1]);

  // Create an object of stencil
  Stencil stencil(filename);

  // read data from the file
  stencil.ReadData();

  // run stencil code for 20 iterations
  stencil.RunStencil(20);

  // output data to another file
  std::string output_name(argv[1]);
  output_name = "output_" + output_name;
  stencil.OutputData(output_name);

  return 0;
}