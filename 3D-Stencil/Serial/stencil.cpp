#include "stencil.hpp"

Stencil::Stencil(std::string fn) : data(NULL)
{
  // store the file name in filename variable
  this->filename = fn;
}

// clean up the memory allocations
Stencil::~Stencil()
{
  // free up data array
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      delete [] data[i][j];
    }
    delete [] data[i];
  }
  delete [] data;
}

// This function will read the data from file name and store it inside the object
// nx, ny, nz are the dimensions of the 3D array
// and the data goes into this->data
Stencil::ReadData()
{
  // local varibles
  int i, j, k;

  // create input file stream
  std::ifstream ifile(filename.c_str());

  // check if the file is open
  if(ifile.is_open()) {
    // the first line should include the dimensions of 3D array
    ifile >> nx >> ny >> nz;

    // check to see if data is already allocated
    // for now let's just exit
    if(data != NULL) {
      std::cerr << "Data is already been allocated!" << std::endl;
      return -1;
    }

    // allocate the 3D array based on the size
    data = new float**[nx];
    for(i=0; i<nx; i++) {
      data[i] = new float*[ny];
      for(j=0; j<ny; j++) {
        data[i][j] = new float[nz];
      }
    }

    // read the data
    for(i=0; i<nx; i++) {
      for(j=0; j<ny; j++) {
        for(k=0; k<nz; k++) {
          ifile >> data[i][j][k];
        }
      }
    }

    // close the file once we are done!
    ifile.close();
  } else {
    std::cerr << "Couldn't open input file" << std::endl;
  }
}