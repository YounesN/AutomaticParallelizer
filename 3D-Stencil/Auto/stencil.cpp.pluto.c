#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))

#include <iostream>
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

  /* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5, t6, t7, t8;
  int lb, ub, lbp, ubp, lb2, ub2;
  register int lbv, ubv;
  /* Start of CLooG code */
  if ((iterations >= 2) && (nx >= 1) && (ny >= 1) && (nz >= 1)) {
    for (t1=0;t1<=floord(2*iterations+nx-4,32);t1++) {
      lbp=std::max(ceild(t1,2),ceild(32*t1-iterations+2,32));
      ubp=std::min(std::min(floord(iterations+nx-2,32),floord(32*t1+nx+31,64)),double(t1));
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
      for (t2=lbp;t2<=ubp;t2++) {
        for (t3=std::max(ceild(32*t2-nx-30,32),double(t1-t2));t3<=std::min(std::min(floord(iterations+ny-2,32),floord(32*t2+ny+30,32)),floord(32*t1-32*t2+ny+31,32));t3++) {
          for (t4=std::max(std::max(ceild(32*t2-nx-30,32),ceild(32*t3-ny-30,32)),double(t1-t2));t4<=std::min(std::min(std::min(floord(iterations+nz-2,32),floord(32*t2+nz+30,32)),floord(32*t3+nz+30,32)),floord(32*t1-32*t2+nz+31,32));t4++) {
            for (t5=std::max(std::max(std::max(32*t1-32*t2,32*t2-nx),32*t3-ny),32*t4-nz);t5<=std::min(std::min(std::min(std::min(iterations-2,32*t2+30),32*t3+30),32*t4+30),32*t1-32*t2+31);t5++) {
              for (t6=std::max(32*t2,t5+1);t6<=std::min(32*t2+31,t5+nx);t6++) {
                for (t7=std::max(32*t3,t5+1);t7<=std::min(32*t3+31,t5+ny);t7++) {
                  lbv=std::max(32*t4,t5+1);
                  ubv=std::min(32*t4+31,t5+nz);
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    A[t5+1][(-t5+t6)][(-t5+t7)][(-t5+t8)] =
                      (A[t5][(-t5+t6)][(-t5+t7)][(-t5+t8)-1] +
                      A[t5][(-t5+t6)][(-t5+t7)][(-t5+t8)+1] +
                      A[t5][(-t5+t6)][(-t5+t7)-1][(-t5+t8)] +
                      A[t5][(-t5+t6)][(-t5+t7)+1][(-t5+t8)] +
                      A[t5][(-t5+t6)-1][(-t5+t7)][(-t5+t8)] +
                      A[t5][(-t5+t6)+1][(-t5+t7)][(-t5+t8)] -
                      6 * A[t5][(-t5+t6)][(-t5+t7)][(-t5+t8)]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
/* End of CLooG code */
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
