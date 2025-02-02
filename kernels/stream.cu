#pragma once
#include "stream.h"
#include "stdio.h"
#include <stdlib.h>
#include <math.h>

#include "deriv.h"

__global__ void check()
{
	printf("\n thread x:%i y:%i, information copied from device:\n", threadIdx.x, threadIdx.y);
	printf("Ra= %f Pr=%f \n", dev.Ra, dev.Pr);
	printf("Rav= %f K=%f Le=%f \n", dev.Rav, dev.K, dev.Le);
	printf("psi1= %f psi2= %f psiS= %f\n", dev.psi1, dev.psi2, dev.psiS);
	printf("Sc1= %f Sc2=%f \n", dev.Sc11, dev.Sc22);
	printf("dim= %i \n", dev.dim);
	printf("hx= %f hy=%f hz=%f \n", dev.hx, dev.hy, dev.hz);
	printf("tau= %20.16f  \n", dev.tau);
	printf("tau_p= %20.16f  \n", dev.tau_p);
	printf("nx= %i ny=%i nz=%i N=%i \n", dev.nx, dev.ny, dev.nz, dev.N);
	printf("Lx= %f Ly=%f Lz=%f \n", dev.Lx, dev.Ly, dev.Lz);
	printf("offset= %i offset2=%i \n", dev.offset, dev.offset2);
	printf("grav_x= %f grav_y= %f \n", dev.grav_x, dev.grav_y);
	printf("vibr_x= %f vibr_y= %f \n", dev.vibr_x, dev.vibr_y);
	printf("density_x= %f density_y= %f \n", dev.density_x, dev.density_y);
	printf("xbc= %i \n", (int)dev.xbc);
	printf("A= %f Ca= %f Gr = %f\n", dev.A, dev.Ca, dev.Gr);

	printf("\n");
}

__global__ void swap_one(double* f_old, double* f_new)
{
	unsigned int l = blockIdx.x * blockDim.x + threadIdx.x;
	if (l < dev.N)
		f_old[l] = f_new[l];
}

__global__ void swap_two(double* f_old, double* f_new, double* f2_old, double* f2_new)
{
	unsigned int l = blockIdx.x * blockDim.x + threadIdx.x;
	if (l < dev.N)
	{
		f_old[l] = f_new[l];
		f2_old[l] = f2_new[l];
	}	
}

__global__ void swap_three(double* f_old, double* f_new, double* f2_old, double* f2_new, double* f3_old, double* f3_new)
{
	unsigned int l = blockIdx.x * blockDim.x + threadIdx.x;
	if (l < dev.N)
	{
		f_old[l] = f_new[l];
		f2_old[l] = f2_new[l];
		f3_old[l] = f3_new[l];
	}
}





namespace stream
{
	#define VX_ dy1(l, ksi)
	#define VY_ -dx1(l, ksi)
	__global__ void vorticity(double* omega_new, double* omega, double* ksi, double* T, double* C)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		auto InnerComputng = [&](unsigned int l)
		{
			return omega[l] + dev.tau * (
				(dx1(l, ksi) * dy1(l, omega) - dy1(l, ksi) * dx1(l, omega)) //nonlinear term
				+ (dx2(l, omega) + dy2(l, omega)) /** dev.Pr*/

				+ dev.grav_y * dev.Ra / dev.Pr * (dx1(l, T) - dev.density_x) 
				- dev.grav_x * dev.Ra / dev.Pr * (dy1(l, T) - dev.density_y)

				+ dev.grav_y * dev.Ra / dev.Pr * dev.K * (dx1(l, C) - dev.density_x)
				- dev.grav_x * dev.Ra / dev.Pr * dev.K * (dy1(l, C) - dev.density_y)
				);            
		};


		if (i <= dev.nx && j <= dev.ny && l < dev.N)
		{
			/*	INNER	*/
			if (i > 0 && i < dev.nx && j > 0 && j < dev.ny)
			{
				omega_new[l] = InnerComputng(l);
			}
			else
			{
				if (j == 0 && (i > 0 && i < dev.nx))
				{
					omega_new[l] = -0.5 / (dev.hy * dev.hy) * (8.0 * ksi[l + dev.offset] - ksi[l + dev.offset * 2]);
					return;
				}
				else if (j == dev.ny && (i > 0 && i < dev.nx))
				{
					omega_new[l] = -0.5 / (dev.hy * dev.hy) * (8.0 * ksi[l - dev.offset] - ksi[l - dev.offset * 2]);
					return;
				}

				if (dev.xbc == closed)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
						omega_new[l] = -0.5 / (dev.hx * dev.hx) * (8.0 * ksi[l + 1] - ksi[l + 2]);
					if (i == dev.nx && (j > 0 && j < dev.ny))
						omega_new[l] = -0.5 / (dev.hx * dev.hx) * (8.0 * ksi[l - 1] - ksi[l - 2]);
					return;
				}
				else if (dev.xbc == periodic)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						unsigned int ll = dev.nx - 1 + dev.offset * j;
						omega_new[l] = InnerComputng(ll);
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						unsigned int ll = 1 + dev.offset * j;
						omega_new[l] = InnerComputng(ll);
						return;
					}
				}
				{
					omega_new[l] = 0;
					omega_new[l] = 0;
				}

			}
		}

	}

	__global__ void make_velocity_from_stream(double* ksi, double* vx, double* vy)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		if (i <= dev.nx && j <= dev.ny && l < dev.N)
		{
			/*	INNER	*/
			if (i > 0 && i < dev.nx && j > 0 && j < dev.ny)
			{
				vx[l] = dy1(l, ksi);
				vy[l] = -dx1(l, ksi);
			}
			else
			{
				if (j == 0 && (i > 0 && i < dev.nx))
				{
					vx[l] = vy[l] = 0.0;
					return;
				}
				else if (j == dev.ny && (i > 0 && i < dev.nx))
				{
					vx[l] = vy[l] = 0.0;
					return;
				}

				if (dev.xbc == closed)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
						vx[l] = vy[l] = 0.0;
					if (i == dev.nx && (j > 0 && j < dev.ny))
						vx[l] = vy[l] = 0.0;
					return;
				}
				else if (dev.xbc == periodic)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						int ll = dev.nx - 1 + dev.offset * j;
						vx[l] = dy1(ll, ksi);
						vy[l] = -dx1(ll, ksi);
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						int ll = 1 + dev.offset * j;
						vx[l] = dy1(ll, ksi);
						vy[l] = -dx1(ll, ksi);
						return;
					}
				}
				{
					vx[l] = 0;
					vy[l] = 0;
				}

			}
		}
	}

	__global__ void poisson_stream(double* ksi_new, double* ksi, double* omega)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		double tau = 0.2 * dev.hx * dev.hy;

		if (i <= dev.nx && j <= dev.ny && l < dev.N)
		{
			/*	INNER	*/
			if (i > 0 && i < dev.nx && j > 0 && j < dev.ny)
			{
				ksi_new[l] = ksi[l] + tau * (dx2(l, ksi) + dy2(l, ksi) + omega[l]);
			}
			else
			{
				if (j == 0 && (i > 0 && i < dev.nx))
				{
					ksi_new[l] = 0.0;
					return;
				}
				else if (j == dev.ny && (i > 0 && i < dev.nx))
				{
					ksi_new[l] = 0.0;
					return;
				}

				if (dev.xbc == closed)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
						ksi_new[l] = 0.0;
						//ksi_new[l] = dx1_eq_0_forward(l, ksi);
					if (i == dev.nx && (j > 0 && j < dev.ny))
						ksi_new[l] = 0.0;
						//ksi_new[l] = dx1_eq_0_back(l, ksi);
					return;
				}
				else if (dev.xbc == periodic)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						int ll = dev.nx - 1 + dev.offset * j;
						ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi));
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						int ll = 1 + dev.offset * j;
						ksi_new[l] = ksi[ll] + tau * (dx2(ll, ksi) + dy2(ll, ksi));
						return;
					}
				}
				{
					ksi_new[l] = 0;
				}

			}
		}

	}

	__global__ void disturb(unsigned int i, unsigned int j, double* f, double value)
	{
		f[i + dev.offset * j] = value;
	}

	__global__ void temperature_2d(double* T, double* T0, double* ksi)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		if (i <= dev.nx && j <= dev.ny && l < dev.N)
		{
			/*	INNER	*/
			if (i > 0 && i < dev.nx && j > 0 && j < dev.ny)
			{
				T[l] = T0[l]
					+ dev.tau * (
						- dy1(l, ksi) * dx1(l, T0) + dx1(l, ksi) * dy1(l, T0)
						+ (VX_ * dev.density_x + VY_ * dev.density_y)

						+ (dx2(l, T0) + dy2(l, T0)) / dev.Pr
						);
				return;
			}

			else
			{
				if (j == 0)
				{
					T[l] = 0.0;
					return;
				}
				else if (j == dev.ny)
				{
					T[l] = 0.0;
					return;
				}

				if (dev.xbc == closed)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						T[l] = dx1_eq_0_forward(l, T0);
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						T[l] = dx1_eq_0_back(l, T0);
						return;
					}

				}
				else if (dev.xbc == periodic)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						int ll = dev.nx - 1 + dev.offset * j;
						T[l] = T0[ll];
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						int ll = 1 + dev.offset * j;
						T[l] = T0[ll];
						return;
					}

				}

				T[l] = 0;
			}
		}
	}

	__global__ void concentration_2d(double* C, double* C0, double* ksi)
	{
		unsigned int i = threadIdx.x + blockIdx.x * blockDim.x;
		unsigned int j = threadIdx.y + blockIdx.y * blockDim.y;
		unsigned int l = i + dev.offset * j;

		if (i <= dev.nx && j <= dev.ny && l < dev.N)
		{
			/*	INNER	*/
			if (i > 0 && i < dev.nx && j > 0 && j < dev.ny)
			{
				C[l] = C0[l]
					+ dev.tau * (
						-dy1(l, ksi) * dx1(l, C0) + dx1(l, ksi) * dy1(l, C0)
						+ (VX_ * dev.density_x + VY_ * dev.density_y)

						+ (dx2(l, C0) + dy2(l, C0)) / (dev.Le * dev.Pr)
						);
				return;
			}

			else
			{
				if (j == 0)
				{
					C[l] = dy1_eq_0_up(l, C0);
					return;
				}
				else if (j == dev.ny)
				{
					C[l] = dy1_eq_0_down(l, C0);
					return;
				}

				if (dev.xbc == closed)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						C[l] = dx1_eq_0_forward(l, C0);
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						C[l] = dx1_eq_0_back(l, C0);
						return;
					}

				}
				else if (dev.xbc == periodic)
				{
					if (i == 0 && (j > 0 && j < dev.ny))
					{
						int ll = dev.nx - 1 + dev.offset * j;
						C[l] = C0[ll];
						return;
					}
					if (i == dev.nx && (j > 0 && j < dev.ny))
					{
						int ll = 1 + dev.offset * j;
						C[l] = C0[ll];
						return;
					}

				}

				C[l] = 0;
			}
		}
	}
}

