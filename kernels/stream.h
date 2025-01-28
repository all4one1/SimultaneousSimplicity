#pragma once
#include "cuda_runtime.h"
#include <device_launch_parameters.h>
#include "../types_project.h"


extern __constant__ Configuration dev;


__global__ void check();
__global__ void swap_one(double* f_old, double* f_new);
__global__ void swap_two(double* f_old, double* f_new, double* f2_old, double* f2_new);
__global__ void swap_three(double* f_old, double* f_new, double* f2_old, double* f2_new, double* f3_old, double* f3_new);

namespace stream
{
	__global__ void vorticity(double* omega_new, double* omega, double* ksi, double* T, double* C);
	__global__ void make_velocity_from_stream(double* ksi, double* vx, double* vy);
	__global__ void temperature_2d(double* T, double* T0, double* ksi);
	__global__ void poisson_stream(double* ksi, double* ksi0, double* omega);
	__global__ void disturb(unsigned int i, unsigned int j, double* f, double val);
	__global__ void concentration_2d(double* C, double* C0, double* ksi);
}