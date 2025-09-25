#pragma once
#include "fv_projection.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

#include "fv_fluxes.h"


namespace fv_projection_cuda
{
	__global__ void temperature_2d(double* T, double* T0, double* ksi)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		double rhs = 0.0;

		if (i <= dev.nx && j <= dev.ny && l < dev.N) // ???
		{
			//x
			if (i > 0)
			{
				rhs += west(T0, l);
			}
			else
			{
				if		(dev.xbc == closed)	  rhs += west(T0, l, 2, 0);
				else if (dev.xbc == periodic) rhs += west(T0, l, 2, 0);
			}

			if (i < dev.nx - 1)
			{
				rhs += east(T0, l);
			}
			else
			{
				if		(dev.xbc == closed)	  rhs += east(T0, l, 2, 0);
				else if (dev.xbc == periodic) rhs += east(T0, l, 2, 0);
			}

			//y
			if (j > 0)
			{
				rhs += south(T0, l);
			}
			else
			{
				rhs += south(T0, l, 2, 0);
			}

			if (j < dev.ny - 1)
			{
				rhs += north(T0, l);
			}
			else
			{
				rhs += north(T0, l, 2, 0);
			}

			T[l] = T0[l] + dev.tau / dev.dV * rhs;
		}
	}
	


}









//__device__ double laplaceX(double* f, unsigned int l, Side s, double bc_val)
//{
//	if (s == Side::center)		return dev.Sx * ((f[l + 1] - f[l]) / dev.hx - ((f[l] - f[l - 1]) / dev.hx));
//	else if (s == Side::west)	return dev.Sx * ((f[l + 1] - f[l]) / dev.hx - ((f[l] - f[l - 1]) / dev.hx));
//	else if (s == Side::east)	return (bc_val);
//}
//
//
//__device__ double laplace(double* f, unsigned int l, unsigned int i, unsigned int j)
//{
//	double lapl = 0.0;
//	lapl += dev.Sx * ((f[l + 1] - f[l]) / dev.hx - ((f[l] - f[l - 1]) / dev.hx));
//	lapl += dev.Sy * ((f[l + dev.offset] - f[l]) / dev.hy - ((f[l] - f[l - dev.offset]) / dev.hy));
//	return lapl;
//}