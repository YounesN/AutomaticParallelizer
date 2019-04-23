#include <iostream>

class Stencil {
  // private members
  std::string filename;
  int nx, ny;
  float **A;
  float **B;
  int padding;
  int iterations;

  // delete a 2D array
  void delete2DArray(float **arr);

  // allocate 2D array with some padding
  void allocate2DArray(float ***arr);
public:
  // constructor of the class
  Stencil(std::string fn, int it);

  // desctructor of the class
  ~Stencil();

  // This function will read the data from file name and store it inside the object
  void ReadData();

  // run stencil for t iterations
  void RunStencil();

  // write data to output file
  void OutputData(std::string output_name);
};