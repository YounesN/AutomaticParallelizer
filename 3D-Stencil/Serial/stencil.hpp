#include <iostream>

class Stencil {
  // private members
  std::string filename;
  int nx, ny, nz;
  float ***data;

public:
  // constructor of the class
  Stencil(std::string fn);

  // desctructor of the class
  ~Stencil();

  // This function will read the data from file name and store it inside the object
  void ReadData();
};