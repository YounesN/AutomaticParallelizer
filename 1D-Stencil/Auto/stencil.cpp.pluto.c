#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

#include "stencil.hpp"
#include <fstream>
#include <cstring>

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
  delete1DArray(A);
}

// delete a 1D array
void Stencil::delete1DArray(float *arr)
{
  // free up data array
  delete [] arr;
}

// allocate a 1D array
void Stencil::allocate1DArray(float **arr)
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
  (*arr) = new float[nx+pad];
  memset((*arr), 0.0, sizeof(float) * (nx+pad));
  
}

// This function will read the data from file name and store it inside the object
// nx, ny, nz are the dimensions of the 3D array
// and the data goes into this->data
void Stencil::ReadData()
{
  // local varibles
  int i;

  // create input file stream
  std::ifstream ifile(filename.c_str());

  // check if the file is open
  if(ifile.is_open()) {
    // the first line should include the dimensions of 3D array
    ifile >> nx;

    // allocate both A0 and ANext arrays
    allocate1DArray(&A);
    allocate1DArray(&B);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      ifile >> A[i];
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

/* Copyright (C) 1991-2018 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */
/* We do not support C11 <threads.h>.  */
  int t1, t2, t3, t4, t5;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((iterations >= 2) && (nx >= 1)) {
  for (t1=0;t1<=floord(3*iterations+nx-5,32);t1++) {
    lbp=max(ceild(2*t1,3),ceild(32*t1-iterations+2,32));
    ubp=min(min(floord(2*iterations+nx-3,32),floord(64*t1+nx+63,96)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      if (t1 <= floord(96*t2-nx-1,64)) {
        if ((nx+1)%2 == 0) {
          A[nx] = B[nx];;
        }
      }
      if (nx == 1) {
        for (t3=16*t2;t3<=min(min(iterations-2,16*t2+14),32*t1-32*t2+31);t3++) {
          B[1] = (A[1 -1] + A[1 +1] + A[1]) * (1.0/3.0);;
          A[1] = B[1];;
        }
      }
      for (t3=max(ceild(32*t2-nx,2),32*t1-32*t2);t3<=min(min(min(floord(32*t2-nx+30,2),iterations-2),16*t2-1),32*t1-32*t2+31);t3++) {
        for (t4=32*t2;t4<=2*t3+nx;t4++) {
          B[(-2*t3+t4)] = (A[(-2*t3+t4)-1] + A[(-2*t3+t4)+1] + A[(-2*t3+t4)]) * (1.0/3.0);;
          A[(-2*t3+t4-1)] = B[(-2*t3+t4-1)];;
        }
        A[nx] = B[nx];;
      }
      for (t3=max(ceild(32*t2-nx+31,2),32*t1-32*t2);t3<=min(min(iterations-2,16*t2-1),32*t1-32*t2+31);t3++) {
        for (t4=32*t2;t4<=32*t2+31;t4++) {
          B[(-2*t3+t4)] = (A[(-2*t3+t4)-1] + A[(-2*t3+t4)+1] + A[(-2*t3+t4)]) * (1.0/3.0);;
          A[(-2*t3+t4-1)] = B[(-2*t3+t4-1)];;
        }
      }
      if (nx >= 2) {
        for (t3=16*t2;t3<=min(min(floord(32*t2-nx+30,2),iterations-2),32*t1-32*t2+31);t3++) {
          B[1] = (A[1 -1] + A[1 +1] + A[1]) * (1.0/3.0);;
          for (t4=2*t3+2;t4<=2*t3+nx;t4++) {
            B[(-2*t3+t4)] = (A[(-2*t3+t4)-1] + A[(-2*t3+t4)+1] + A[(-2*t3+t4)]) * (1.0/3.0);;
            A[(-2*t3+t4-1)] = B[(-2*t3+t4-1)];;
          }
          A[nx] = B[nx];;
        }
      }
      for (t3=max(ceild(32*t2-nx+31,2),16*t2);t3<=min(min(iterations-2,16*t2+14),32*t1-32*t2+31);t3++) {
        B[1] = (A[1 -1] + A[1 +1] + A[1]) * (1.0/3.0);;
        for (t4=2*t3+2;t4<=32*t2+31;t4++) {
          B[(-2*t3+t4)] = (A[(-2*t3+t4)-1] + A[(-2*t3+t4)+1] + A[(-2*t3+t4)]) * (1.0/3.0);;
          A[(-2*t3+t4-1)] = B[(-2*t3+t4-1)];;
        }
      }
      if ((t1 >= ceild(3*t2-1,2)) && (t2 <= floord(iterations-17,16))) {
        B[1] = (A[1 -1] + A[1 +1] + A[1]) * (1.0/3.0);;
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

  ofile << nx << std::endl;
  for(int i=padding ; i<nx+padding; i++) {
    ofile << A[i] << " ";
  }
  ofile << std::endl;
}
