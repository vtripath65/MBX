
/******************************************************************************
Copyright 2019 The Regents of the University of California.
All Rights Reserved.

Permission to copy, modify and distribute any part of this Software for
educational, research and non-profit purposes, without fee, and without
a written agreement is hereby granted, provided that the above copyright
notice, this paragraph and the following three paragraphs appear in all
copies.

Those desiring to incorporate this Software into commercial products or
use for commercial purposes should contact the:
Office of Innovation & Commercialization
University of California, San Diego
9500 Gilman Drive, Mail Code 0910
La Jolla, CA 92093-0910
Ph: (858) 534-5815
FAX: (858) 534-7345
E-MAIL: invent@ucsd.edu

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY
OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS. THE UNIVERSITY OF CALIFORNIA MAKES NO
REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND, EITHER IMPLIED OR
EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT THE USE OF THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR OTHER RIGHTS.
******************************************************************************/

#include "poly_2b_A1B2Z2_A1B2Z2_deg4_vh2o_revPBE0_def2svpd.h"
#include "poly_2b_A1B2Z2_A1B2Z2_deg4_grad_vh2o_revPBE0_def2svpd_cu.h"
#include <cuda.h>
#include <iostream>
/**
 * @file poly_2b_A1B2Z2_A1B2Z2_deg4_nograd_vh2o_revPBE0_def2svpd.cpp
 * @brief Contains the implementation of the polynomials without gradients for symmetry A1B2Z2_A1B2Z2
 */

/**
 * @namespace mbnrg_A1B2Z2_A1B2Z2_deg4
 * @brief Encloses the structure of the polynomial for symmetry A1B2Z2_A1B2Z2
 */

#define safeMalloc(pntr, size) \
do { \
    cudaError_t err = cudaMalloc(pntr, size); \
    if (err != cudaSuccess) { \
        std::cout << "CUDA malloc failed: " << cudaGetErrorString(err) << std::endl; \
   } \
} while(0)

struct device_memory_grad{
  double* devx;
  double* deva;
  double* devg;
};

device_memory_grad devmemgrad;


//__constant__ device_memory_grad d_devmemgrad;

namespace mbnrg_A1B2Z2_A1B2Z2_deg4 {


double poly_A1B2Z2_A1B2Z2_deg4_vh2o_revPBE0_def2svpd::eval_direct(const double x[31], const double a[1208], double g[31])
{

    double energy = 0.0;
    double *denergy;

    std::fill(g, g + 31, 0.0);

//    std::cout << "Launching the kernel" << std::endl;

    safeMalloc((void**)&devmemgrad.devx, 31*sizeof(double));
    safeMalloc((void**)&devmemgrad.deva, 1208*sizeof(double));
    safeMalloc((void**)&devmemgrad.devg, 31*sizeof(double));
    safeMalloc((void**)&denergy, sizeof(double));

    cudaMemcpy(devmemgrad.devx, x, 31*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devmemgrad.deva, a, 1208*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devmemgrad.devg, g, 31*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(denergy, &energy, sizeof(double), cudaMemcpyHostToDevice);

//    cudaMemcpyToSymbol(d_devmemgrad.devx, devmemgrad.devx, sizeof(double*));
//    cudaMemcpyToSymbol(d_devmemgrad.deva, devmemgrad.deva, sizeof(double*));
//    cudaMemcpyToSymbol(d_devmemgrad.devg, devmemgrad.devg, sizeof(double*));

    kernel_2B<<<14,256>>>(denergy, devmemgrad.devx, devmemgrad.deva, devmemgrad.devg);
//    cudaDeviceSynchronize();

    cudaMemcpy(&energy, denergy, sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(g, devmemgrad.devg, 31*sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(devmemgrad.devx);
    cudaFree(devmemgrad.deva);
    cudaFree(devmemgrad.devg);
    cudaFree(denergy);

//    cudaDeviceReset();

//    cudaFree(d_devmemgrad.devx);
//    cudaFree(d_devmemgrad.deva);
//    cudaFree(d_devmemgrad.devg);

//   std::cout << "enegy_cuda = " << energy << std::endl;

    return energy;
};

}; // namespace mbnrg_A1B2Z2_A1B2Z2_deg4
