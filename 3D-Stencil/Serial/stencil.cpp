#include "stencil.hpp"
#include <fstream>

Stencil::Stencil(std::string fn) : A0(NULL), ANext(NULL)
{
  // store the file name in filename variable
  this->filename = fn;
    
  // we are going to add some padding to simplify stencil operations
  padding = 1;
}

// clean up the memory allocations
Stencil::~Stencil()
{
  // free both arrays
  delete3DArray(A0);
  delete3DArray(ANext);
}

// delete a 3D array
void Stencil::delete3DArray(float ***arr)
{
  // free up data array
  for(int i=0; i<nx; i++) {
    for(int j=0; j<ny; j++) {
      delete [] arr[i][j];
    }
    delete [] arr[i];
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
  (*arr) = new float**[nx + pad];
  for(int i=0; i<nx+pad; i++) {
    (*arr)[i] = new float*[ny + pad];
    for(int j=0; j<ny+pad; j++) {
      (*arr)[i][j] = new float[nz + pad];

      // set every elemento to zero
      memset((*arr)[i][j], 0.0, sizeof(float) * (nz+pad));
    }
  }
}

void Stencil::swap3DArray(float ***arr1, float ***arr2)
{
  float ***tmp = arr1;
  arr1 = arr2;
  arr2 = tmp;
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
    allocate3DArray(&A0);
    allocate3DArray(&ANext);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        for(k=padding; k<nz+padding; k++) {
          ifile >> A0[i][j][k];
        }
      }
    }

    // close the file once we are done!
    ifile.close();
  } else {
    std::cerr << "Couldn't open input file" << std::endl;
  }
}

void Stencil::RunStencil(int t)
{
  int i, j, k, it;

  for(it=0; it<t; it++) {
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        for(k=padding; k<nz+padding; k++) {
          ANext[i][j][k] = (A0[i][j][k-1] + A0[i][j][k+1] + 
                            A0[i][j-1][k] + A0[i][j+1][k] +
                            A0[i-1][j][k] + A0[i+1][j][k] -
                            6 * A0[i][j][k]);
        }
      }
    }
    swap3DArray(A0, ANext);
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
        ofile << A0[i][j][k] << " ";
      }
      ofile << std::endl;
    }
    ofile << std::endl;
  }
}