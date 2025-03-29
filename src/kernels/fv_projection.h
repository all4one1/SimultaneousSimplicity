#pragma once
#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include "../types_project.h"


extern __constant__ Configuration dev;


__global__ void check();
__global__ void swap_one(double* f_old, double* f_new);
__global__ void swap_two(double* f_old, double* f_new, double* f2_old, double* f2_new);
__global__ void swap_three(double* f_old, double* f_new, double* f2_old, double* f2_new, double* f3_old, double* f3_new);


namespace fv_projection_cuda
{


}

