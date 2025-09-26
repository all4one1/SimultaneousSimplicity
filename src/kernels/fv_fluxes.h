#pragma once

//
//
//__device__ double west(double* f, unsigned int l, int bc = 0, double bc_val = 0)
//{
//	if (bc == 0)		return (f[l] - f[l - 1]) / dev.hx;
//	else if (bc == 1)	return (f[l] - bc_val) / (0.5 * dev.hx);
//	else if (bc == 2)	return (bc_val);
//	else if (bc == 3)   return (f[l] - f[l + (dev.nx - 1)]) / dev.hx;
//}
//
//__device__ double east(double* f, unsigned int l, int bc = 0, double bc_val = 0)
//{
//	if (bc == 0)		return (f[l + 1] - f[l]) / dev.hx;
//	else if (bc == 1)	return (bc_val - f[l]) / (0.5 * dev.hx);
//	else if (bc == 2)	return (bc_val);
//	else if (bc == 3)   return (f[(dev.nx - 1)] - f[l]) / dev.hx;
//}
//
//__device__ double south(double* f, unsigned int l, int bc = 0, double bc_val = 0)
//{
//	if (bc == 0)		return (f[l] - f[l - dev.offset]) / dev.hy;
//	else if (bc == 1)	return (f[l] - bc_val) / (0.5 * dev.hy);
//	else if (bc == 2)	return (bc_val);
//	else if (bc == 3)   return (f[l] - f[l + dev.offset * (dev.ny - 1)]) / dev.hy;
//}
//
//__device__ double north(double* f, unsigned int l, int bc = 0, double bc_val = 0)
//{
//	if (bc == 0)		return (f[l + dev.offset] - f[l]) / dev.hy;
//	else if (bc == 1)	return (bc_val - f[l]) / (0.5 * dev.hy);
//	else if (bc == 2)	return (bc_val);
//	else if (bc == 3)   return (dev.offset * f[(dev.ny - 1)] - f[l]) / dev.hy;
//}
