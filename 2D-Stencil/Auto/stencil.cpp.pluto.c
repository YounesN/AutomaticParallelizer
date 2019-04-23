#include <omp.h>
#include <math.h>
#include "stencil.hpp"
#include <fstream>
#include <cstring>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))


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
  int t1, t2, t3, t4, t5, t6;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((iterations >= 2) && (nx >= 1) && (ny >= 1)) {
  for (t1=ceild(padding-31,32);t1<=floord(3*iterations+padding+nx-6,32);t1++) {
    lbp=max(ceild(32*t1-iterations+2,32),ceild(64*t1+padding-31,96));
    ubp=min(min(floord(2*iterations+padding+nx-4,32),floord(64*t1+padding+nx+62,96)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(ceild(32*t2-nx-30,32),ceild(64*t1-64*t2+padding-31,32));t3<=min(min(floord(32*t2+ny+30,32),floord(2*iterations+padding+ny-4,32)),floord(64*t1-64*t2+padding+ny+62,32));t3++) {
        if ((t1 <= floord(96*t2-padding-nx,64)) && (t2 >= ceild(32*t3+nx-ny+1,32))) {
          if ((padding+nx)%2 == 0) {
            for (t6=max(32*t3,32*t2-nx+1);t6<=min(32*t3+31,32*t2-nx+ny);t6++) {
              A[(padding+nx-1)][(-32*t2+t6+padding+nx-1)] = B[(padding+nx-1)][(-32*t2+t6+padding+nx-1)];;
            }
          }
        }
        if ((t1 <= floord(64*t2+32*t3-padding-ny,64)) && (t2 <= floord(32*t3+nx-ny,32))) {
          if ((padding+ny)%2 == 0) {
            for (t5=max(32*t2,32*t3-ny+1);t5<=min(32*t2+31,32*t3+nx-ny);t5++) {
              A[(-32*t3+t5+padding+ny-1)][(padding+ny-1)] = B[(-32*t3+t5+padding+ny-1)][(padding+ny-1)];;
            }
          }
        }
        if ((nx >= 2) && (ny == 1) && (t2 == t3)) {
          for (t4=max(ceild(32*t2-padding,2),32*t1-32*t2);t4<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t4++) {
            B[padding][padding] = (A[padding][padding-1] + A[padding][padding+1] + A[padding-1][padding] + A[padding+1][padding] + A[padding][padding]) * (1.0/5.0);;
            for (t5=2*t4+padding+1;t5<=2*t4+padding+nx-1;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][padding] = B[(-2*t4+t5-1)][padding];;
            }
            A[(padding+nx-1)][padding] = B[(padding+nx-1)][padding];;
          }
        }
        if ((ny == 1) && (t2 == t3)) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t4<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t4++) {
            B[padding][padding] = (A[padding][padding-1] + A[padding][padding+1] + A[padding-1][padding] + A[padding+1][padding] + A[padding][padding]) * (1.0/5.0);;
            for (t5=2*t4+padding+1;t5<=32*t2+31;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][padding] = B[(-2*t4+t5-1)][padding];;
            }
          }
        }
        if (nx == 1) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t4<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=max(32*t3,2*t4+padding);t6<=min(32*t3+31,2*t4+padding+ny-1);t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t6=max(32*t3,2*t4+padding+1);t6<=min(32*t3+31,2*t4+padding+ny);t6++) {
              A[padding][(-2*t4+t6-1)] = B[padding][(-2*t4+t6-1)];;
            }
          }
        }
        if (ny == 1) {
          for (t4=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=2*t4+padding+nx-1;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][padding] = B[(-2*t4+t5-1)][padding];;
            }
            A[(padding+nx-1)][padding] = B[(padding+nx-1)][padding];;
          }
        }
        if (ny == 1) {
          for (t4=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][padding] = B[(-2*t4+t5-1)][padding];;
            }
          }
        }
        for (t4=max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t4<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=2*t4+padding+nx-1;t5++) {
            for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
            A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
          }
          for (t6=32*t3;t6<=2*t4+padding+ny;t6++) {
            A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
          }
        }
        for (t4=max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=32*t2+31;t5++) {
            for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
            A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
          }
        }
        for (t4=max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=2*t4+padding+nx-1;t5++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
          }
          for (t6=32*t3;t6<=32*t3+31;t6++) {
            A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
          }
        }
        for (t4=max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=32*t2+31;t5++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if ((nx >= 2) && (ny >= 2) && (t2 == t3)) {
          for (t4=max(ceild(32*t2-padding,2),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=2*t4+padding;t6<=2*t4+padding+ny-1;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=2*t4+padding+nx-1;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=2*t4+padding+ny-1;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
              A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
            }
            for (t6=2*t4+padding+1;t6<=2*t4+padding+ny;t6++) {
              A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if ((ny >= 2) && (t2 == t3)) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t4<=min(min(floord(32*t2-padding-ny+31,2),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=2*t4+padding;t6<=2*t4+padding+ny-1;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=32*t2+31;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=2*t4+padding+ny-1;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
              A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
            }
          }
        }
        if (nx >= 2) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=2*t4+padding+nx-1;t5++) {
              for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
              A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
            }
            for (t6=32*t3;t6<=2*t4+padding+ny;t6++) {
              A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
            }
          }
        }
        for (t4=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
            B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
          }
          for (t5=2*t4+padding+1;t5<=32*t2+31;t5++) {
            for (t6=32*t3;t6<=2*t4+padding+ny-1;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
            A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
          }
        }
        if ((nx >= 2) && (t2 == t3)) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=2*t4+padding;t6<=32*t2+31;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=2*t4+padding+nx-1;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=32*t2+31;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
            }
            for (t6=2*t4+padding+1;t6<=32*t2+31;t6++) {
              A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if (t2 == t3) {
          for (t4=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=2*t4+padding;t6<=32*t2+31;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=32*t2+31;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=32*t2+31;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
            }
          }
        }
        if (nx >= 2) {
          for (t4=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
            }
            for (t5=2*t4+padding+1;t5<=2*t4+padding+nx-1;t5++) {
              for (t6=32*t3;t6<=32*t3+31;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
            }
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
            }
          }
        }
        for (t4=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t6=32*t3;t6<=32*t3+31;t6++) {
            B[padding][(-2*t4+t6)] = (A[padding][(-2*t4+t6)-1] + A[padding][(-2*t4+t6)+1] + A[padding-1][(-2*t4+t6)] + A[padding+1][(-2*t4+t6)] + A[padding][(-2*t4+t6)]) * (1.0/5.0);;
          }
          for (t5=2*t4+padding+1;t5<=32*t2+31;t5++) {
            for (t6=32*t3;t6<=32*t3+31;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if (ny >= 2) {
          for (t4=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=2*t4+padding+nx-1;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=2*t4+padding+ny-1;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
              A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
            }
            for (t6=2*t4+padding+1;t6<=2*t4+padding+ny;t6++) {
              A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if (ny >= 2) {
          for (t4=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
            for (t5=32*t2;t5<=32*t2+31;t5++) {
              B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
              for (t6=2*t4+padding+1;t6<=2*t4+padding+ny-1;t6++) {
                B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
                A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
              }
              A[(-2*t4+t5-1)][(padding+ny-1)] = B[(-2*t4+t5-1)][(padding+ny-1)];;
            }
          }
        }
        for (t4=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=2*t4+padding+nx-1;t5++) {
            B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
            for (t6=2*t4+padding+1;t6<=32*t3+31;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
          }
          for (t6=2*t4+padding+1;t6<=32*t3+31;t6++) {
            A[(padding+nx-1)][(-2*t4+t6-1)] = B[(padding+nx-1)][(-2*t4+t6-1)];;
          }
        }
        for (t4=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t4<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),iterations-2),32*t1-32*t2+31);t4++) {
          for (t5=32*t2;t5<=32*t2+31;t5++) {
            B[(-2*t4+t5)][padding] = (A[(-2*t4+t5)][padding-1] + A[(-2*t4+t5)][padding+1] + A[(-2*t4+t5)-1][padding] + A[(-2*t4+t5)+1][padding] + A[(-2*t4+t5)][padding]) * (1.0/5.0);;
            for (t6=2*t4+padding+1;t6<=32*t3+31;t6++) {
              B[(-2*t4+t5)][(-2*t4+t6)] = (A[(-2*t4+t5)][(-2*t4+t6)-1] + A[(-2*t4+t5)][(-2*t4+t6)+1] + A[(-2*t4+t5)-1][(-2*t4+t6)] + A[(-2*t4+t5)+1][(-2*t4+t6)] + A[(-2*t4+t5)][(-2*t4+t6)]) * (1.0/5.0);;
              A[(-2*t4+t5-1)][(-2*t4+t6-1)] = B[(-2*t4+t5-1)][(-2*t4+t6-1)];;
            }
          }
        }
        if ((t1 >= ceild(96*t2-padding-31,64)) && (t2 <= min(floord(2*iterations+padding-35,32),t3-1))) {
          if ((padding+1)%2 == 0) {
            for (t6=32*t3;t6<=min(32*t3+31,32*t2+ny+30);t6++) {
              B[padding][(-32*t2+t6+padding-31)] = (A[padding][(-32*t2+t6+padding-31)-1] + A[padding][(-32*t2+t6+padding-31)+1] + A[padding-1][(-32*t2+t6+padding-31)] + A[padding+1][(-32*t2+t6+padding-31)] + A[padding][(-32*t2+t6+padding-31)]) * (1.0/5.0);;
            }
          }
        }
        if ((t1 >= ceild(64*t2+32*t3-padding-31,64)) && (t2 >= t3) && (t3 <= floord(2*iterations+padding-35,32))) {
          if ((padding+1)%2 == 0) {
            for (t5=max(32*t2,32*t3+31);t5<=min(32*t2+31,32*t3+nx+30);t5++) {
              B[(-32*t3+t5+padding-31)][padding] = (A[(-32*t3+t5+padding-31)][padding-1] + A[(-32*t3+t5+padding-31)][padding+1] + A[(-32*t3+t5+padding-31)-1][padding] + A[(-32*t3+t5+padding-31)+1][padding] + A[(-32*t3+t5+padding-31)][padding]) * (1.0/5.0);;
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

  ofile << nx << " " << ny << std::endl;
  for(int i=padding ; i<nx+padding; i++) {
    for(int j=padding; j<ny+padding; j++) {
      ofile << A[i][j] << " ";
    }
    ofile << std::endl;
  }
}
