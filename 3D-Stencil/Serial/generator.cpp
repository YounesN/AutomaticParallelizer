#include <iostream>
#include <stdlib.h>     /* atoi */
#include <cstdlib>      /* rand and srand */
#include <ctime>        /* time */

int main(int argc, char *argv[])
{
  // check the number of parameters
  if(argc < 5) {
    // Tell the user how to run the program
    std::cerr << "Usage: " << argv[0] << " <filename> <nx> <ny> <nz>" << std::endl;
    return -1;
  }

  // output object
  std::ofstream ofile(argv[1]);

  // check to see if the file is open
  if(!ofile.is_open()) {
    std::cerr << "Output file is not open!" << std::endl;
    return -1;
  }

  // dimension sizes
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  int nz = atoi(argv[4]);
  float max_number = 100;

  // provide seed number
  srand (static_cast <unsigned> (time(0)));

  // output the first line which is the dimensions
  ofile << nx << ny << nz << std::endl;
  for(int i=0 ; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      for(int k=0; k<nz; k++) {
        float r = static_cast <float> (rand()) /
          (static_cast <float> (RAND_MAX/max_number));
        ofile << r << " ";
      }
      ofile << std::endl;
    }
    ofile << std::endl;
  }

  // close the file
  ofile.close();

  return 0;
}