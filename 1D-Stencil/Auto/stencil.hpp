#include <iostream>

class Stencil {
  // private members
  std::string filename;
  int nx;
  float *A;
  float *B;
  int padding;
  int iterations;

  // delete a 1D array
  void delete1DArray(float *arr);

  // allocate 1D array with some padding
  void allocate1DArray(float **arr);
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