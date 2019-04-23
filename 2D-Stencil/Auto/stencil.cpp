#include "stencil.hpp"
#include <fstream>
#include <cstring>
#include <omp.h>

Stencil::Stencil(std::string fn, int it)
{
  // store inputs into member variables
  this->filename = fn;
  this->iterations = it;
    
  // we are going to add some padding to simplify stencil operations
  padding = 1;

  // Set A and B to NULL
  A = NULL;
  B = NULL;
}

// clean up the memory allocations
Stencil::~Stencil()
{
  // free array
  delete2DArray(A);
}

// delete a 2D array
void Stencil::delete2DArray(float **arr)
{
  // free up data array
  for(int i=0; i<nx; i++) {
    delete [] arr[i];
  }
  delete [] arr;
}

// allocate a 2D array
void Stencil::allocate2DArray(float ***arr)
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
  (*arr) = new float*[nx+pad];
  for(int i=0; i<nx+pad; i++) {
    (*arr)[i] = new float[ny + pad];
      // set every elemento to zero
      memset((*arr)[i], 0.0, sizeof(float) * (ny+pad));
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
    allocate2DArray(&A);
    allocate2DArray(&B);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        ifile >> A[i][j];
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

#pragma scop
  for(it=0; it<iterations-1; it++) {
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
          B[i][j] = (A[i][j-1] + A[i][j+1] +
                     A[i-1][j] + A[i+1][j] +
                     A[i][j]) * (1.0/5.0);
      }
    }
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
          A[i][j] = B[i][j];
      }
    }
  }
#pragma endscop
}

void Stencil::OutputData(std::string output_name)
{
  std::ofstream ofile;
  ofile.open(output_name.c_str());

  ofile << nx << " " << ny << std::endl;
  for(int i=padding ; i<nx+padding; i++) {
    for(int j=padding; j<ny+padding; j++) {
      ofile << A[i][j] << " ";
    }
    ofile << std::endl;
  }
}