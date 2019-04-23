#include <iostream>
#include <stdlib.h>     /* atoi */
#include <cstdlib>      /* rand and srand */
#include <ctime>        /* time */
#include <fstream>      /* ofstream */

int main(int argc, char *argv[])
{
  // check the number of parameters
  if(argc < 3) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <filename> <dim>" << std::endl;
    return -1;
  }

  // output object
  std::ofstream ofile;
  ofile.open(argv[1]);

  // check to see if the file is open
  if(!ofile.is_open()) {
    std::cerr << "Output file is not open!" << std::endl;
    return -1;
  }

  // dimension sizes
  int nx = atoi(argv[2]);
  float max_number = 100;

  // provide seed number
  srand (static_cast <unsigned> (time(0)));

  // output the first line which is the dimensions
  ofile << nx << std::endl;
  for(int i=0 ; i<nx; i++) {
    float r = static_cast <float> (rand()) /
      (static_cast <float> (RAND_MAX/max_number));
    ofile << r << " ";
  }
  ofile << std::endl;

  // close the file
  ofile.close();

  return 0;
}
