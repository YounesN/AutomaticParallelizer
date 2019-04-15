#include <iostream>

class Stencil {
  // private members
  std::string filename;
  int nx, ny, nz;
  float ***A0;
  float ***ANext;
  int padding;

  // delete a 3D array
  void delete3DArray(float ***arr);

  // allocate 3D array with some padding
  void allocate3DArray(float ****arr);

  // swap the content of arrays
  void swap3DArray(float ***arr1, float ***arr2);
public:
  // constructor of the class
  Stencil(std::string fn);

  // desctructor of the class
  ~Stencil();

  // This function will read the data from file name and store it inside the object
  void ReadData();

  // run stencil for t iterations
  void RunStencil(int t);

  // write data to output file
  void OutputData(std::string output_name);
};