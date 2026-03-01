#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "types_project.h"
#include <iostream>
#include <vector>
#include <array>
using std::cout;
using std::endl;

struct Boundary
{
	double value = 0.0;
	MathBoundary type = MathBoundary::Neumann;
	double* ptr = nullptr;
	Boundary() {};
	Boundary(MathBoundary t, double v = 0.0) : type(t), value(v) {};

	double& operator()(unsigned int q) { return ptr[q]; }
};

//using walls = std::array<Boundary, 6>;
struct DomainBounds
{
	std::array<Boundary, 6> boundaries;
	int dim = 1;
	double* v = nullptr;
	int nx, ny, nz, offset, offset2;

	DomainBounds(MathBoundary def_type = MathBoundary::Neumann, double def_value = 0.0,
		double hx = 0, double hy = 0, double hz = 0)
	{
		(*this)[Side::west] = Boundary(def_type, def_value);
		(*this)[Side::east] = Boundary(def_type, def_value);
		(*this)[Side::south] = Boundary(def_type, def_value);
		(*this)[Side::north] = Boundary(def_type, def_value);
		(*this)[Side::front] = Boundary(def_type, def_value);
		(*this)[Side::back] = Boundary(def_type, def_value);
		nx = ny = nz = offset = offset2 = 0;
	}
	Boundary& operator[](Side side) {
		return boundaries[static_cast<size_t>(side)];
	}
	const Boundary& operator[](Side side) const {
		return boundaries[static_cast<size_t>(side)];
	}
};

struct Variable
{
	double* f = nullptr;
	DomainBounds bounds;
	Configuration conf;

	Variable()
	{

	};
	Variable(double* ptr, Configuration conf_) : f(ptr), conf(conf_)
	{

	}
	double& operator[](unsigned int i) { return f[i]; }
	operator double* () const { return f; }
};

template<Side side>
__device__ double value(Variable* f_, unsigned int l, bool wall = false)
{
	const auto& f = *f_;
	const auto& hx = f_->conf.hx;
	const auto& hy = f_->conf.hy;
	const auto& hz = f_->conf.hz;
					 
	const auto& nx = f_->conf.nx;
	const auto& ny = f_->conf.ny;
	const auto& nz = f_->conf.nz;

	const auto& offset = f_->conf.offset;
	const auto& offset2 = f_->conf.offset2;
	const auto& border_value = f_->bounds[side].value;
	const auto& type = f_->bounds[side].type;


	if constexpr (side == Side::west)
	{
		if (!wall)
			return (f[l] + f[l - 1]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] - 0.5 * hx * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l] + f[l + (nx - 1)]) * 0.5;
	}
	else if constexpr (side == Side::east)
	{
		if (!wall)
			return (f[l + 1] + f[l]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] + 0.5 * hx * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l - (nx - 1)] + f[l]) * 0.5;
	}
	else if constexpr (side == Side::south)
	{
		if (!wall)
			return (f[l] + f[l - offset]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] - 0.5 * hy * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l] + f[l + offset * (ny - 1)]) * 0.5;
	}
	else if constexpr (side == Side::north) {
		if (!wall)
			return (f[l + offset] + f[l]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] + 0.5 * hy * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l - offset * (ny - 1)] + f[l]) * 0.5;
	}
	else if constexpr (side == Side::front)
	{
		if (!wall)
			return (f[l] + f[l - offset2]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] - 0.5 * hz * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l] + f[l + offset2 * (nz - 1)]) * 0.5;
	}
	else if constexpr (side == Side::back) {
		if (!wall)
			return (f[l + offset2] + f[l]) * 0.5;
		else if (type == MathBoundary::Dirichlet)	  return (border_value);
		else if (type == MathBoundary::Neumann)	  return f[l] + 0.5 * hz * border_value;
		else if (type == MathBoundary::periodic)	  return (f[l - offset2 * (nz - 1)] + f[l]) * 0.5;
	}
	return 0;
}

template<Side side>
__device__ double deriv(Variable* f_, unsigned int l, bool wall = false)
{
	const auto& f = *f_;
	const auto& hx = f_->conf.hx;
	const auto& hy = f_->conf.hy;
	const auto& hz = f_->conf.hz;

	const auto& nx = f_->conf.nx;
	const auto& ny = f_->conf.ny;
	const auto& nz = f_->conf.nz;

	const auto& offset = f_->conf.offset;
	const auto& offset2 = f_->conf.offset2;
	const auto& border_value = f_->bounds[side].value;
	const auto& type = f_->bounds[side].type;

	if constexpr (side == Side::west) {
		if (!wall)
			return (f[l] - f[l - 1]) / hx;
		else if (type == MathBoundary::Dirichlet)    return (f[l] - border_value) / (0.5 * hx);
		else if (type == MathBoundary::Neumann)      return (border_value);
		else if (type == MathBoundary::periodic)     return (f[l] - f[l + (nx - 1)]) / hx;
	}
	else if constexpr (side == Side::east) {
		if (!wall)
			return (f[l + 1] - f[l]) / hx;
		else if (type == MathBoundary::Dirichlet)   return (border_value - f[l]) / (0.5 * hx);
		else if (type == MathBoundary::Neumann)     return (border_value);
		else if (type == MathBoundary::periodic)    return (f[l - (nx - 1)] - f[l]) / hx;
	}
	else if constexpr (side == Side::south)
	{
		if (!wall)
			return (f[l] - f[l - offset]) / hy;
		else if (type == MathBoundary::Dirichlet)   return (f[l] - border_value) / (0.5 * hy);
		else if (type == MathBoundary::Neumann)     return (border_value);
		else if (type == MathBoundary::periodic)    return (f[l] - f[l + offset * (ny - 1)]) / hy;
	}
	else if constexpr (side == Side::north) {
		if (!wall)
			return (f[l + offset] - f[l]) / hy;
		else if (type == MathBoundary::Dirichlet)   return (border_value - f[l]) / (0.5 * hy);
		else if (type == MathBoundary::Neumann)     return (border_value);
		else if (type == MathBoundary::periodic)    return (f[l - offset * (ny - 1)] - f[l]) / hy;
	}
	else if constexpr (side == Side::front)
	{
		if (!wall)
			return (f[l] - f[l - offset]) / hy;
		else if (type == MathBoundary::Dirichlet)   return (f[l] - border_value) / (0.5 * hz);
		else if (type == MathBoundary::Neumann)     return (border_value);
		else if (type == MathBoundary::periodic)    return (f[l] - f[l + offset2 * (ny - 1)]) / hz;
	}
	else if constexpr (side == Side::back) {
		if (!wall)
			return (f[l + offset] - f[l]) / hy;
		else if (type == MathBoundary::Dirichlet)   return (border_value - f[l]) / (0.5 * hz);
		else if (type == MathBoundary::Neumann)     return (border_value);
		else if (type == MathBoundary::periodic)    return (f[l - offset2 * (ny - 1)] - f[l]) / hz;
	}
	return 0; // дл€ подавлени€ предупреждений компил€тора
}


void init_bc(std::vector<Variable> &all, Configuration &config)
{
	if (config.domain == closed_box)
	{
		//all[field_T][Side::west] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_T])[east] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_T])[south] = { MathBoundary::Dirichlet, 1.0 };
		//(*all[field_T])[north] = { MathBoundary::Dirichlet, 0.0 };

		//(*all[field_vx])[west] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx])[east] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx])[south] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx])[north] = { MathBoundary::Dirichlet, 0.0 };

		//(*all[field_vy])[west] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy])[east] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy])[south] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy])[north] = { MathBoundary::Dirichlet, 0.0 };

		//(*all[field_p])[west] = { MathBoundary::complex, 0.0 };
		//(*all[field_p])[east] = { MathBoundary::complex, 0.0 };
		//(*all[field_p])[south] = { MathBoundary::complex, 0.0 };
		//(*all[field_p])[north] = { MathBoundary::complex, 0.0 };

		//(*all[field_p_prime])[west] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_p_prime])[east] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_p_prime])[south] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_p_prime])[north] = { MathBoundary::Neumann, 0.0 };

		//(*all[field_vx_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vx_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

		//(*all[field_vy_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
		//(*all[field_vy_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

		//(*all[field_C])[west] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_C])[east] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_C])[south] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_C])[north] = { MathBoundary::Neumann, 0.0 };

		//(*all[field_mu])[west] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_mu])[east] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_mu])[south] = { MathBoundary::Neumann, 0.0 };
		//(*all[field_mu])[north] = { MathBoundary::Neumann, 0.0 };
	}
}



__constant__ Configuration config;
//
//
//int main()
//{
//
//
//	Variable* host;
//	Variable* device;
//
//	host = new Variable;
//	cudaMalloc((void**)&device, sizeof(Variable));
//
//
//	cudaMemcpy(host, device, sizeof(Variable), cudaMemcpyHostToDevice);
//
//
//
//
//	return 0;
//}
