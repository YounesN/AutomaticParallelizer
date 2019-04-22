#include "stencil.hpp"
#include <fstream>
#include <cstring>

Stencil::Stencil(std::string fn, int it) : A(NULL)
{
  // store inputs into member variables
  this->filename = fn;
  this->iterations = it;
    
  // we are going to add some padding to simplify stencil operations
  padding = 1;
}

// clean up the memory allocations
Stencil::~Stencil()
{
  // free array
  delete3DArray(A);
}

// delete a 3D array
void Stencil::delete3DArray(float ***arr)
{
  // free up data array
  for(int it=0; it<iterations; it++) {
    for(int i=0; i<nx; i++) {
      delete [] arr[it][i];
    }
    delete [] arr[it];
  }
  delete [] arr;
}

// allocate a 3D array
void Stencil::allocate3DArray(float ****arr)
{
  // check to see if data is already allocated
  // for now let's just exit
  if((*arr) != NULL) {
    std::cerr << "Data is already been allocated!" << std::endl;
    return;
  }

  // add padding on each side
  // so multiply the padding by 2
  int pad = padding * 2;

  // allocate the 3D array based on the size
  // and initialize to 0.0
  (*arr) = new float**[iterations];
  for(int it=0; it<iterations; it++) {
    (*arr)[it] = new float*[nx + pad];
    for(int i=0; i<nx+pad; i++) {
      (*arr)[it][i] = new float[ny + pad];
      
      memset((*arr)[it][i], 0.0, sizeof(float) * (ny + pad));
    }
  }
}

// This function will read the data from file name and store it inside the object
// nx, ny, nz are the dimensions of the 3D array
// and the data goes into this->data
void Stencil::ReadData()
{
  // local varibles
  int i, j;

  // create input file stream
  std::ifstream ifile(filename.c_str());

  // check if the file is open
  if(ifile.is_open()) {
    // the first line should include the dimensions of 3D array
    ifile >> nx >> ny;

    // allocate both A0 and ANext arrays
    allocate3DArray(&A);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
	ifile >> A[0][i][j];
      }
    }

    // close the file once we are done!
    ifile.close();
  } else {
    std::cerr << "Couldn't open input file" << std::endl;
  }
}

void Stencil::RunStencil()
{
  int i, j, it;
  for(it=0; it<iterations-1; it++) {
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
	      A[it+1][i][j] = (A[it][i][j-1] + A[it][i][j+1] +
          A[it][i-1][j] + A[it][i+1][j] - 4 * A[it][i][j]);
      }
    }
  }
}

void Stencil::OutputData(std::string output_name)
{
  std::ofstream ofile;
  ofile.open(output_name.c_str());
  
  ofile << nx << " " << ny << " " << std::endl;
  for(int i=padding ; i<nx+padding; i++) {
    for(int j=padding; j<ny+padding; j++) {
      ofile << A[iterations-1][i][j] << " ";
    }
    ofile << std::endl;
  }
}
