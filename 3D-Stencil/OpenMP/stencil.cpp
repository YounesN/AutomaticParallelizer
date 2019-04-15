#include "stencil.hpp"
#include <fstream>
#include <cstring>
#include <omp.h>

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
  delete4DArray(A);
}

// delete a 3D array
void Stencil::delete4DArray(float ****arr)
{
  // free up data array
  for(int it=0; it<iterations; it++) {
    for(int i=0; i<nx; i++) {
      for(int j=0; j<ny; j++) {
        delete [] arr[it][i][j];
      }
      delete [] arr[it][i];
    }
    delete [] arr[it];
  }
  delete [] arr;
}

// allocate a 3D array
void Stencil::allocate4DArray(float *****arr)
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
  (*arr) = new float***[iterations];
  for(int it=0; it<iterations; it++) {
    (*arr)[it] = new float**[nx + pad];
    for(int i=0; i<nx+pad; i++) {
      (*arr)[it][i] = new float*[ny + pad];
      for(int j=0; j<ny+pad; j++) {
        (*arr)[it][i][j] = new float[nz + pad];

        // set every elemento to zero
        memset((*arr)[it][i][j], 0.0, sizeof(float) * (nz+pad));
      }
    }
  }
}

// This function will read the data from file name and store it inside the object
// nx, ny, nz are the dimensions of the 3D array
// and the data goes into this->data
void Stencil::ReadData()
{
  // local varibles
  int i, j, k;

  // create input file stream
  std::ifstream ifile(filename.c_str());

  // check if the file is open
  if(ifile.is_open()) {
    // the first line should include the dimensions of 3D array
    ifile >> nx >> ny >> nz;

    // allocate both A0 and ANext arrays
    allocate4DArray(&A);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        for(k=padding; k<nz+padding; k++) {
          ifile >> A[0][i][j][k];
        }
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
  int i, j, k, it;

#ifdef _OPENMP
    #pragma omp parallel for default(shared) private(it, i, j, k)
#endif
  for(it=0; it<iterations-1; it++) {
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        for(k=padding; k<nz+padding; k++) {
          A[it+1][i][j][k] = (A[it][i][j][k-1] + A[it][i][j][k+1] + 
                            A[it][i][j-1][k] + A[it][i][j+1][k] +
                            A[it][i-1][j][k] + A[it][i+1][j][k] -
                            6 * A[it][i][j][k]);
        }
      }
    }
  }
}

void Stencil::OutputData(std::string output_name)
{
  std::ofstream ofile;
  ofile.open(output_name.c_str());

  ofile << nx << " " << ny << " " << nz << std::endl;
  for(int i=padding ; i<nx+padding; i++) {
    for(int j=padding; j<ny+padding; j++) {
      for(int k=padding; k<nz+padding; k++) {
        ofile << A[iterations-1][i][j][k] << " ";
      }
      ofile << std::endl;
    }
    ofile << std::endl;
  }
}