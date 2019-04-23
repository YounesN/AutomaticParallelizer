#include <omp.h>
#include <math.h>
#include "stencil.hpp"
#include <fstream>
#include <cstring>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    (std::max(double(x), double(y)))
#define min(x,y)    (std::min(double(x), double(y)))


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
  (*arr) = new float**[nx+pad];
  for(int i=0; i<nx+pad; i++) {
    (*arr)[i] = new float*[ny + pad];
    for(int j=0; j<ny+pad; j++) {
      (*arr)[i][j] = new float[nz + pad];

      // set every elemento to zero
      memset((*arr)[i][j], 0.0, sizeof(float) * (nz+pad));
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
    allocate3DArray(&A);
    allocate3DArray(&B);

    // read the data
    for(i=padding; i<nx+padding; i++) {
      for(j=padding; j<ny+padding; j++) {
        for(k=padding; k<nz+padding; k++) {
          ifile >> A[i][j][k];
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
  int t1, t2, t3, t4, t5, t6, t7, t8;
 int lb, ub, lbp, ubp, lb2, ub2;
 register int lbv, ubv;
/* Start of CLooG code */
if ((iterations >= 2) && (nx >= 1) && (ny >= 1) && (nz >= 1)) {
  for (t1=ceild(padding-31,32);t1<=floord(3*iterations+padding+nx-6,32);t1++) {
    lbp=max(ceild(32*t1-iterations+2,32),ceild(64*t1+padding-31,96));
    ubp=min(min(floord(2*iterations+padding+nx-4,32),floord(64*t1+padding+nx+62,96)),t1);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=max(ceild(32*t2-nx-30,32),ceild(64*t1-64*t2+padding-31,32));t3<=min(min(floord(32*t2+ny+30,32),floord(2*iterations+padding+ny-4,32)),floord(64*t1-64*t2+padding+ny+62,32));t3++) {
        for (t4=max(max(ceild(32*t2-nx-30,32),ceild(32*t3-ny-30,32)),ceild(64*t1-64*t2+padding-31,32));t4<=min(min(min(floord(32*t2+nz+30,32),floord(32*t3+nz+30,32)),floord(2*iterations+padding+nz-4,32)),floord(64*t1-64*t2+padding+nz+62,32));t4++) {
          if ((t1 <= floord(96*t2-padding-nx,64)) && (t2 >= max(ceild(32*t3+nx-ny+1,32),ceild(32*t4+nx-nz+1,32)))) {
            if ((padding+nx)%2 == 0) {
              for (t7=max(32*t3,32*t2-nx+1);t7<=min(32*t3+31,32*t2-nx+ny);t7++) {
                lbv=max(32*t4,32*t2-nx+1);
                ubv=min(32*t4+31,32*t2-nx+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-32*t2+t7+padding+nx-1)][(-32*t2+t8+padding+nx-1)] = B[(padding+nx-1)][(-32*t2+t7+padding+nx-1)][(-32*t2+t8+padding+nx-1)];;
                }
              }
            }
          }
          if ((t1 <= floord(64*t2+32*t3-padding-ny,64)) && (t2 <= floord(32*t3+nx-ny,32)) && (t3 >= ceild(32*t4+ny-nz+1,32))) {
            if ((padding+ny)%2 == 0) {
              for (t6=max(32*t2,32*t3-ny+1);t6<=min(32*t2+31,32*t3+nx-ny);t6++) {
                lbv=max(32*t4,32*t3-ny+1);
                ubv=min(32*t4+31,32*t3-ny+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-32*t3+t6+padding+ny-1)][(padding+ny-1)][(-32*t3+t8+padding+ny-1)] = B[(-32*t3+t6+padding+ny-1)][(padding+ny-1)][(-32*t3+t8+padding+ny-1)];;
                }
              }
            }
          }
          if ((t1 <= floord(64*t2+32*t4-padding-nz,64)) && (t2 <= floord(32*t4+nx-nz,32)) && (t3 <= floord(32*t4+ny-nz,32))) {
            if ((padding+nz)%2 == 0) {
              for (t6=max(32*t2,32*t4-nz+1);t6<=min(32*t2+31,32*t4+nx-nz);t6++) {
                for (t7=max(32*t3,32*t4-nz+1);t7<=min(32*t3+31,32*t4+ny-nz);t7++) {
                  A[(-32*t4+t6+padding+nz-1)][(-32*t4+t7+padding+nz-1)][(padding+nz-1)] = B[(-32*t4+t6+padding+nz-1)][(-32*t4+t7+padding+nz-1)][(padding+nz-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (ny >= 2) && (nz == 1) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(ceild(32*t2-padding,2),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((ny >= 2) && (nz == 1) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-ny+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
            }
          }
          if ((nx >= 2) && (nz == 1) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((nz == 1) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
            }
          }
          if ((nx >= 2) && (nz == 1) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((nz == 1) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
            }
          }
          if ((nx >= 2) && (nz == 1) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((nz == 1) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                B[padding][(-2*t5+t7)][padding] = (A[padding][(-2*t5+t7)][padding-1] + A[padding][(-2*t5+t7)][padding+1] + A[padding][(-2*t5+t7)-1][padding] + A[padding][(-2*t5+t7)+1][padding] + A[padding-1][(-2*t5+t7)][padding] + A[padding+1][(-2*t5+t7)][padding] + A[padding][(-2*t5+t7)][padding]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
            }
          }
          if ((nx >= 2) && (ny == 1) && (t2 == t3)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              lbv=max(32*t4,2*t5+padding);
              ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][padding][(-2*t5+t8)] = (A[padding][padding][(-2*t5+t8)-1] + A[padding][padding][(-2*t5+t8)+1] + A[padding][padding-1][(-2*t5+t8)] + A[padding][padding+1][(-2*t5+t8)] + A[padding-1][padding][(-2*t5+t8)] + A[padding+1][padding][(-2*t5+t8)] + A[padding][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=max(32*t4,2*t5+padding);
                ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                lbv=max(32*t4,2*t5+padding+1);
                ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][padding][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][padding][(-2*t5+t8-1)];;
                }
              }
              lbv=max(32*t4,2*t5+padding+1);
              ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][padding][(-2*t5+t8-1)] = B[(padding+nx-1)][padding][(-2*t5+t8-1)];;
              }
            }
          }
          if ((ny == 1) && (t2 == t3)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t5++) {
              lbv=max(32*t4,2*t5+padding);
              ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][padding][(-2*t5+t8)] = (A[padding][padding][(-2*t5+t8)-1] + A[padding][padding][(-2*t5+t8)+1] + A[padding][padding-1][(-2*t5+t8)] + A[padding][padding+1][(-2*t5+t8)] + A[padding-1][padding][(-2*t5+t8)] + A[padding+1][padding][(-2*t5+t8)] + A[padding][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=max(32*t4,2*t5+padding);
                ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                lbv=max(32*t4,2*t5+padding+1);
                ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][padding][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][padding][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nx == 1) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=max(32*t3,2*t5+padding);t7<=min(32*t3+31,2*t5+padding+ny-1);t7++) {
                lbv=max(32*t4,2*t5+padding);
                ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t7=max(32*t3,2*t5+padding+1);t7<=min(32*t3+31,2*t5+padding+ny);t7++) {
                lbv=max(32*t4,2*t5+padding+1);
                ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[padding][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[padding][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (nz == 1) && (t3 == t4)) {
            for (t5=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((ny >= 2) && (nz == 1) && (t3 == t4)) {
            for (t5=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
            }
          }
          if ((nz == 1) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if ((nz == 1) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                B[(-2*t5+t6)][padding][padding] = (A[(-2*t5+t6)][padding][padding-1] + A[(-2*t5+t6)][padding][padding+1] + A[(-2*t5+t6)][padding-1][padding] + A[(-2*t5+t6)][padding+1][padding] + A[(-2*t5+t6)-1][padding][padding] + A[(-2*t5+t6)+1][padding][padding] + A[(-2*t5+t6)][padding][padding]) * (1.0/7.0);;
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
            }
          }
          if (ny == 1) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=max(32*t4,2*t5+padding);
                ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                lbv=max(32*t4,2*t5+padding+1);
                ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][padding][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][padding][(-2*t5+t8-1)];;
                }
              }
              lbv=max(32*t4,2*t5+padding+1);
              ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][padding][(-2*t5+t8-1)] = B[(padding+nx-1)][padding][(-2*t5+t8-1)];;
              }
            }
          }
          if (ny == 1) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=max(32*t4,2*t5+padding);
                ubv=min(32*t4+31,2*t5+padding+nz-1);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                lbv=max(32*t4,2*t5+padding+1);
                ubv=min(32*t4+31,2*t5+padding+nz);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][padding][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][padding][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nz == 1) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if (nz == 1) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
                A[(-2*t5+t6-1)][(padding+ny-1)][padding] = B[(-2*t5+t6-1)][(padding+ny-1)][padding];;
              }
            }
          }
          if (nz == 1) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][padding] = B[(padding+nx-1)][(-2*t5+t7-1)][padding];;
              }
            }
          }
          if (nz == 1) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][padding] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][padding];;
                }
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
            for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
            for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+1,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(ceild(32*t2-padding-nx+32,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (ny >= 2) && (nz >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(ceild(32*t2-padding,2),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-ny+31,2)),floord(32*t2-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (nz >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-ny+31,2),floord(32*t2-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (nz >= 2) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t2-padding-nz+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nz >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nz+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          if ((nx >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nz+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t2 == t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-nz+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          if ((nx >= 2) && (ny >= 2) && (t2 == t3)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t2-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (t2 == t3)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          if ((nx >= 2) && (ny >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-ny+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t2 == t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (t3 >= t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 >= t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (ny >= 2) && (t2 == t3) && (t2 <= t4-1)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (t2 == t3) && (t2 <= t4-1)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-ny+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (t3 <= t4-1)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 <= t4-1) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (nz >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t2-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nz >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-nz+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (nz >= 2) && (t2 == t4)) {
            for (t5=max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t2-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nz >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-nz+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-nz+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
          }
          if ((nx >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nz+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-nx+31,2),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t2 == t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-nz+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nx >= 2) && (t2 == t3)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t2 == t3) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
            }
          }
          if (nx >= 2) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
              }
            }
            for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
          }
          if ((nx >= 2) && (t2 == t3) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),ceild(32*t2-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((t2 == t3) && (t2 == t4)) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),ceild(32*t2-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if ((nx >= 2) && (t2 == t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t2 == t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding;
                ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t2+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if ((nx >= 2) && (t3 >= t4)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t4-padding-1,2),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 >= t4) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding+30,2),floord(32*t4-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if ((nx >= 2) && (t2 == t3) && (t2 <= t4-1)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding-nx+31,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((t2 == t3) && (t2 <= t4-1)) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t2-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(floord(32*t2-padding+30,2),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=2*t5+padding;t7<=32*t2+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t2+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if ((nx >= 2) && (t3 <= t4-1)) {
            for (t5=max(max(max(ceild(32*t2-padding,2),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t3-padding-1,2),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 <= t4-1) {
            for (t5=max(max(max(max(ceild(32*t2-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding+30,2),floord(32*t3-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-2*t5+t7)][(-2*t5+t8)] = (A[padding][(-2*t5+t7)][(-2*t5+t8)-1] + A[padding][(-2*t5+t7)][(-2*t5+t8)+1] + A[padding][(-2*t5+t7)-1][(-2*t5+t8)] + A[padding][(-2*t5+t7)+1][(-2*t5+t8)] + A[padding-1][(-2*t5+t7)][(-2*t5+t8)] + A[padding+1][(-2*t5+t7)][(-2*t5+t8)] + A[padding][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                }
              }
              for (t6=2*t5+padding+1;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if ((ny >= 2) && (nz >= 2) && (t3 == t4)) {
            for (t5=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (nz >= 2) && (t3 == t4)) {
            for (t5=max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-ny+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nz >= 2) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((nz >= 2) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
            }
          }
          if (ny >= 2) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (ny >= 2) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=32*t4;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
            for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t4-padding-1,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              lbv=32*t4;
              ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
                A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
              }
            }
          }
          if ((ny >= 2) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((ny >= 2) && (t3 == t4)) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=2*t5+padding+1;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 == t4) {
            for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=2*t5+padding;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (t3 == t4) {
            for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t3-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=2*t5+padding;
                ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=32*t3+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
              }
            }
          }
          if (ny >= 2) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=2*t5+padding+1;t7<=2*t5+padding+ny;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (ny >= 2) {
            for (t5=max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t4-padding-1,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
                }
                for (t7=2*t5+padding+1;t7<=2*t5+padding+ny-1;t7++) {
                  lbv=32*t4;
                  ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                }
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t4-padding-1,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
            for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t3-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding+30,2)),floord(32*t4-padding-1,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              lbv=32*t4;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                B[(-2*t5+t6)][padding][(-2*t5+t8)] = (A[(-2*t5+t6)][padding][(-2*t5+t8)-1] + A[(-2*t5+t6)][padding][(-2*t5+t8)+1] + A[(-2*t5+t6)][padding-1][(-2*t5+t8)] + A[(-2*t5+t6)][padding+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][padding][(-2*t5+t8)] + A[(-2*t5+t6)+1][padding][(-2*t5+t8)] + A[(-2*t5+t6)][padding][(-2*t5+t8)]) * (1.0/7.0);;
              }
              for (t7=2*t5+padding+1;t7<=32*t3+31;t7++) {
                lbv=32*t4;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nz >= 2) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
              for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nz >= 2) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t3-padding-ny+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nz >= 2) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t2-padding-nx+31,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                lbv=2*t5+padding+1;
                ubv=2*t5+padding+nz;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if (nz >= 2) {
            for (t5=max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding-nz+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
              for (t6=32*t2;t6<=32*t2+31;t6++) {
                for (t7=32*t3;t7<=32*t3+31;t7++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                  lbv=2*t5+padding+1;
                  ubv=2*t5+padding+nz-1;
#pragma ivdep
#pragma vector always
                  for (t8=lbv;t8<=ubv;t8++) {
                    B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                    A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                  }
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(padding+nz-1)];;
                }
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t2-padding-nx+31,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                lbv=2*t5+padding+1;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
              lbv=2*t5+padding+1;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
            for (t7=32*t3;t7<=2*t5+padding+ny;t7++) {
              lbv=2*t5+padding+1;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+1,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t3-padding-ny+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=2*t5+padding+ny-1;t7++) {
                B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                lbv=2*t5+padding+1;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
              lbv=2*t5+padding+1;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(padding+ny-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+1,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),floord(32*t2-padding-nx+31,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=2*t5+padding+nx-1;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                lbv=2*t5+padding+1;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
            for (t7=32*t3;t7<=32*t3+31;t7++) {
              lbv=2*t5+padding+1;
              ubv=32*t4+31;
#pragma ivdep
#pragma vector always
              for (t8=lbv;t8<=ubv;t8++) {
                A[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(padding+nx-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
              }
            }
          }
          for (t5=max(max(max(max(ceild(32*t4-padding,2),ceild(32*t2-padding-nx+32,2)),ceild(32*t3-padding-ny+32,2)),ceild(32*t4-padding-nz+32,2)),32*t1-32*t2);t5<=min(min(min(min(floord(32*t2-padding-1,2),floord(32*t3-padding-1,2)),floord(32*t4-padding+30,2)),iterations-2),32*t1-32*t2+31);t5++) {
            for (t6=32*t2;t6<=32*t2+31;t6++) {
              for (t7=32*t3;t7<=32*t3+31;t7++) {
                B[(-2*t5+t6)][(-2*t5+t7)][padding] = (A[(-2*t5+t6)][(-2*t5+t7)][padding-1] + A[(-2*t5+t6)][(-2*t5+t7)][padding+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][padding] + A[(-2*t5+t6)][(-2*t5+t7)+1][padding] + A[(-2*t5+t6)-1][(-2*t5+t7)][padding] + A[(-2*t5+t6)+1][(-2*t5+t7)][padding] + A[(-2*t5+t6)][(-2*t5+t7)][padding]) * (1.0/7.0);;
                lbv=2*t5+padding+1;
                ubv=32*t4+31;
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)] = (A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)-1] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)+1] + A[(-2*t5+t6)][(-2*t5+t7)-1][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)+1][(-2*t5+t8)] + A[(-2*t5+t6)-1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)+1][(-2*t5+t7)][(-2*t5+t8)] + A[(-2*t5+t6)][(-2*t5+t7)][(-2*t5+t8)]) * (1.0/7.0);;
                  A[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)] = B[(-2*t5+t6-1)][(-2*t5+t7-1)][(-2*t5+t8-1)];;
                }
              }
            }
          }
          if ((t1 >= ceild(96*t2-padding-31,64)) && (t2 <= min(min(floord(2*iterations+padding-35,32),t3-1),t4-1))) {
            if ((padding+1)%2 == 0) {
              for (t7=32*t3;t7<=min(32*t3+31,32*t2+ny+30);t7++) {
                lbv=32*t4;
                ubv=min(32*t4+31,32*t2+nz+30);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[padding][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)] = (A[padding][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)-1] + A[padding][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)+1] + A[padding][(-32*t2+t7+padding-31)-1][(-32*t2+t8+padding-31)] + A[padding][(-32*t2+t7+padding-31)+1][(-32*t2+t8+padding-31)] + A[padding-1][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)] + A[padding+1][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)] + A[padding][(-32*t2+t7+padding-31)][(-32*t2+t8+padding-31)]) * (1.0/7.0);;
                }
              }
            }
          }
          if ((t1 >= ceild(64*t2+32*t3-padding-31,64)) && (t2 >= t3) && (t3 <= min(floord(2*iterations+padding-35,32),t4-1))) {
            if ((padding+1)%2 == 0) {
              for (t6=max(32*t2,32*t3+31);t6<=min(32*t2+31,32*t3+nx+30);t6++) {
                lbv=32*t4;
                ubv=min(32*t4+31,32*t3+nz+30);
#pragma ivdep
#pragma vector always
                for (t8=lbv;t8<=ubv;t8++) {
                  B[(-32*t3+t6+padding-31)][padding][(-32*t3+t8+padding-31)] = (A[(-32*t3+t6+padding-31)][padding][(-32*t3+t8+padding-31)-1] + A[(-32*t3+t6+padding-31)][padding][(-32*t3+t8+padding-31)+1] + A[(-32*t3+t6+padding-31)][padding-1][(-32*t3+t8+padding-31)] + A[(-32*t3+t6+padding-31)][padding+1][(-32*t3+t8+padding-31)] + A[(-32*t3+t6+padding-31)-1][padding][(-32*t3+t8+padding-31)] + A[(-32*t3+t6+padding-31)+1][padding][(-32*t3+t8+padding-31)] + A[(-32*t3+t6+padding-31)][padding][(-32*t3+t8+padding-31)]) * (1.0/7.0);;
                }
              }
            }
          }
          if ((t1 >= ceild(64*t2+32*t4-padding-31,64)) && (t2 >= t4) && (t3 >= t4) && (t4 <= floord(2*iterations+padding-35,32))) {
            if ((padding+1)%2 == 0) {
              for (t6=max(32*t2,32*t4+31);t6<=min(32*t2+31,32*t4+nx+30);t6++) {
                for (t7=max(32*t3,32*t4+31);t7<=min(32*t3+31,32*t4+ny+30);t7++) {
                  B[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)][padding] = (A[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)][padding-1] + A[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)][padding+1] + A[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)-1][padding] + A[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)+1][padding] + A[(-32*t4+t6+padding-31)-1][(-32*t4+t7+padding-31)][padding] + A[(-32*t4+t6+padding-31)+1][(-32*t4+t7+padding-31)][padding] + A[(-32*t4+t6+padding-31)][(-32*t4+t7+padding-31)][padding]) * (1.0/7.0);;
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
        ofile << A[i][j][k] << " ";
      }
      ofile << std::endl;
    }
    ofile << std::endl;
  }
}
