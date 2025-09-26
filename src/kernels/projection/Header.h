#pragma once
#include "FromOuterSparse/SparseMatrix.h"
#include "types_project.h"

extern Configuration config;


namespace fv2
{
	struct Boundary
	{
		double value = 0.0;
		MathBoundary type = MathBoundary::Neumann;
		double* v = nullptr;
		Boundary() {};
		Boundary(MathBoundary t, double v = 0.0) : type(t), value(v) {};
		Boundary(MathBoundary t, bool new_array = false, unsigned int n = 0, double fixed_val = NAN) : type(t)
		{
			if (new_array) v = new double[n];
			if (fixed_val == fixed_val)
			{
				for (unsigned int l = 0; l < n; l++)
					v[l] = fixed_val;
			}
		};
		//Boundary(MathBoundary t, bool new_array, unsigned int n, double (*func)(double)) : type(t)
		double& operator()(unsigned int q) { return v[q]; }
	};
	using walls = std::vector<Boundary>;
	struct Variable
	{
		double* f = nullptr;
		walls wall;

		Variable()
		{
			wall.resize(4);
		};
		Variable(double* ptr) : f(ptr)
		{
			wall.resize(4);
		}
		double& operator[](int i) { return f[i]; }
		operator double* () const { return f; }

		void set_walls(Boundary west_, Boundary east_, Boundary south_, Boundary north_)
		{
			wall[west] = west_;
			wall[east] = east_;
			wall[south] = south_;
			wall[north] = north_;
		}
	};

	void init_bc(std::vector<walls*> allwalls)
	{
		int size = 4;
		for (auto& it : allwalls)
			it->resize(size);

		if (config.domain == closed_box)
		{
			(*allwalls[field_T])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_T])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_T])[south] = { MathBoundary::Dirichlet, 1.0 };
			(*allwalls[field_T])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vx])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_p])[west] = { MathBoundary::complex, 0.0 };
			(*allwalls[field_p])[east] = { MathBoundary::complex, 0.0 };
			(*allwalls[field_p])[south] = { MathBoundary::complex, 0.0 };
			(*allwalls[field_p])[north] = { MathBoundary::complex, 0.0 };

			(*allwalls[field_p_prime])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_p_prime])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_p_prime])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_p_prime])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_vx_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_C])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_mu])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[north] = { MathBoundary::Neumann, 0.0 };
		}



		if (config.domain == periodic_cell)
		{
			(*allwalls[field_T])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_T])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_T])[south] = { MathBoundary::Dirichlet, 1.0 };
			(*allwalls[field_T])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vx])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vx])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vx])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vy])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vy])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_p])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_p])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_p])[south] = { MathBoundary::complex, 0.0 };
			(*allwalls[field_p])[north] = { MathBoundary::complex, 0.0 };

			(*allwalls[field_p_prime])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_p_prime])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_p_prime])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_p_prime])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_vx_prime])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vx_prime])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vx_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy_prime])[west] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vy_prime])[east] = { MathBoundary::periodic, 0.0 };
			(*allwalls[field_vy_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[north] = { MathBoundary::Dirichlet, 0.0 };
		}

		if (config.domain == open_tube)
		{
			(*allwalls[field_T])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_T])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_T])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_T])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_vx])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_vx])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_vx])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_p])[west] = { MathBoundary::Dirichlet, config.p_in };
			(*allwalls[field_p])[east] = { MathBoundary::Dirichlet, config.p_out };
			(*allwalls[field_p])[south] = { MathBoundary::complex, 0.0 };
			(*allwalls[field_p])[north] = { MathBoundary::complex, 0.0 };

			(*allwalls[field_p_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_p_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_p_prime])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_p_prime])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_vx_prime])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_vx_prime])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_vx_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vx_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_vy_prime])[west] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[east] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[south] = { MathBoundary::Dirichlet, 0.0 };
			(*allwalls[field_vy_prime])[north] = { MathBoundary::Dirichlet, 0.0 };

			(*allwalls[field_C])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_C])[north] = { MathBoundary::Neumann, 0.0 };

			(*allwalls[field_mu])[west] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[east] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[south] = { MathBoundary::Neumann, 0.0 };
			(*allwalls[field_mu])[north] = { MathBoundary::Neumann, 0.0 };
		}



	}

	template<Side side>
	double deriv(double* f, unsigned int l, bool wall = false, MathBoundary type = MathBoundary::not_boundary, double value = 0)
	{
		if constexpr (side == Side::west) {
			if (!wall)
				return (f[l] - f[l - 1]) / config.hx;
			else
			{
				if (type == MathBoundary::Dirichlet)    return (f[l] - value) / (0.5 * config.hx);
				if (type == MathBoundary::Neumann)      return (value);
				if (type == MathBoundary::periodic)     return (f[l] - f[l + (config.nx - 1)]) / config.hx;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::east) {
			if (!wall)
				return (f[l + 1] - f[l]) / config.hx;
			else
			{
				if (type == MathBoundary::Dirichlet)   return (value - f[l]) / (0.5 * config.hx);
				if (type == MathBoundary::Neumann)     return (value);
				if (type == MathBoundary::periodic)    return (f[l - (config.nx - 1)] - f[l]) / config.hx;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::south)
		{
			if (!wall)
				return (f[l] - f[l - config.offset]) / config.hy;
			else
			{
				if (type == MathBoundary::Dirichlet)   return (f[l] - value) / (0.5 * config.hy);
				if (type == MathBoundary::Neumann)     return (value);
				if (type == MathBoundary::periodic)    return (f[l] - f[l + config.offset * (config.ny - 1)]) / config.hy;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::north) {
			if (!wall)
				return (f[l + config.offset] - f[l]) / config.hy;
			else
			{
				if (type == MathBoundary::Dirichlet)   return (value - f[l]) / (0.5 * config.hy);
				if (type == MathBoundary::Neumann)     return (value);
				if (type == MathBoundary::periodic)    return (f[l - config.offset * (config.ny - 1)] - f[l]) / config.hy;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else {
			static_assert(false, "Invalid Side template parameter");
		}
		return 0; // дл€ подавлени€ предупреждений компил€тора
	}
	template<Side side>
	double value(double* f, unsigned int l, bool wall = false, MathBoundary type = MathBoundary::not_boundary, double value = 0)
	{
		if constexpr (side == Side::west)
		{
			if (!wall)
				return (f[l] + f[l - 1]) * 0.5;
			else
			{
				if (type == MathBoundary::Dirichlet)	 return (value);
				if (type == MathBoundary::Neumann)		 return f[l] - 0.5 * config.hx * value;
				if (type == MathBoundary::periodic)		 return (f[l] + f[l + (config.nx - 1)]) * 0.5;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::east) {
			if (!wall)
				return (f[l + 1] + f[l]) * 0.5;
			else
			{
				if (type == MathBoundary::Dirichlet)	 return (value);
				if (type == MathBoundary::Neumann)		 return f[l] + 0.5 * config.hx * value;
				if (type == MathBoundary::periodic)		 return (f[l - (config.nx - 1)] + f[l]) * 0.5;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::south) {
			if (!wall)
				return (f[l] + f[l - config.offset]) * 0.5;
			else
			{
				if (type == MathBoundary::Dirichlet)	 return (value);
				if (type == MathBoundary::Neumann)		 return f[l] - 0.5 * config.hy * value;
				if (type == MathBoundary::periodic)		 return (f[l] + f[l + config.offset * (config.ny - 1)]) * 0.5;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else if constexpr (side == Side::north) {
			if (!wall)
				return (f[l + config.offset] + f[l]) * 0.5;
			else
			{
				if (type == MathBoundary::Dirichlet)	  return (value);
				if (type == MathBoundary::Neumann)		  return f[l] + 0.5 * config.hy * value;
				if (type == MathBoundary::periodic)		  return (f[l - config.offset * (config.ny - 1)] + f[l]) * 0.5;
				if (type == MathBoundary::array_value) { pause return NAN; };
			}
		}
		else {
			static_assert(false, "Invalid Side template parameter");
		}
		return 0;
	}

	template<Side side>
	double value2(double* f, unsigned int l, bool wall = false, const Boundary& b = Boundary())
	{
		if constexpr (side == Side::west)
		{
			if (!wall)
				return (f[l] + f[l - 1]) * 0.5;
			else if (b.type == MathBoundary::Dirichlet)	  return (b.value);
			else if (b.type == MathBoundary::Neumann)	  return f[l] - 0.5 * config.hx * b.value;
			else if (b.type == MathBoundary::periodic)	  return (f[l] + f[l + (config.nx - 1)]) * 0.5;
			else if (b.type == MathBoundary::array_value) return b.v[l / config.offset];
		}
		else if constexpr (side == Side::east)
		{
			if (!wall)
				return (f[l + 1] + f[l]) * 0.5;
			else if (b.type == MathBoundary::Dirichlet)	  return (b.value);
			else if (b.type == MathBoundary::Neumann)	  return f[l] + 0.5 * config.hx * b.value;
			else if (b.type == MathBoundary::periodic)	  return (f[l - (config.nx - 1)] + f[l]) * 0.5;
			else if (b.type == MathBoundary::array_value) return b.v[l / config.offset];
		}
		else if constexpr (side == Side::south)
		{
			if (!wall)
				return (f[l] + f[l - config.offset]) * 0.5;
			else if (b.type == MathBoundary::Dirichlet)	  return (b.value);
			else if (b.type == MathBoundary::Neumann)	  return f[l] - 0.5 * config.hy * b.value;
			else if (b.type == MathBoundary::periodic)	  return (f[l] + f[l + config.offset * (config.ny - 1)]) * 0.5;
			else if (b.type == MathBoundary::array_value) return b.v[l % config.offset];
		}
		else if constexpr (side == Side::north) {
			if (!wall)
				return (f[l + config.offset] + f[l]) * 0.5;
			else if (b.type == MathBoundary::Dirichlet)	  return (b.value);
			else if (b.type == MathBoundary::Neumann)	  return f[l] + 0.5 * config.hy * b.value;
			else if (b.type == MathBoundary::periodic)	  return (f[l - config.offset * (config.ny - 1)] + f[l]) * 0.5;
			else if (b.type == MathBoundary::array_value) return b.v[l % config.offset];

		}
		else {
			static_assert(false, "Invalid Side template parameter");
			pause
		}
		return 0;
	}
	template<Side side>
	double deriv2(double* f, unsigned int l, bool wall = false, const Boundary& b = Boundary())
	{
		if constexpr (side == Side::west) {
			if (!wall)
				return (f[l] - f[l - 1]) / config.hx;
			else if (b.type == MathBoundary::Dirichlet)    return (f[l] - b.value) / (0.5 * config.hx);
			else if (b.type == MathBoundary::Neumann)      return (b.value);
			else if (b.type == MathBoundary::periodic)     return (f[l] - f[l + (config.nx - 1)]) / config.hx;
			else if (b.type == MathBoundary::array_value)  return (f[l] - b.v[l / config.offset]) / (0.5 * config.hx);
		}
		else if constexpr (side == Side::east) {
			if (!wall)
				return (f[l + 1] - f[l]) / config.hx;
			else if (b.type == MathBoundary::Dirichlet)   return (b.value - f[l]) / (0.5 * config.hx);
			else if (b.type == MathBoundary::Neumann)     return (b.value);
			else if (b.type == MathBoundary::periodic)    return (f[l - (config.nx - 1)] - f[l]) / config.hx;
			else if (b.type == MathBoundary::array_value) return (b.v[l / config.offset] - f[l]) / (0.5 * config.hx);
		}
		else if constexpr (side == Side::south)
		{
			if (!wall)
				return (f[l] - f[l - config.offset]) / config.hy;
			else if (b.type == MathBoundary::Dirichlet)   return (f[l] - b.value) / (0.5 * config.hy);
			else if (b.type == MathBoundary::Neumann)     return (b.value);
			else if (b.type == MathBoundary::periodic)    return (f[l] - f[l + config.offset * (config.ny - 1)]) / config.hy;
			else if (b.type == MathBoundary::array_value) return (f[l] - b.v[l % config.offset]) / (0.5 * config.hy);
		}
		else if constexpr (side == Side::north) {
			if (!wall)
				return (f[l + config.offset] - f[l]) / config.hy;
			else if (b.type == MathBoundary::Dirichlet)   return (b.value - f[l]) / (0.5 * config.hy);
			else if (b.type == MathBoundary::Neumann)     return (b.value);
			else if (b.type == MathBoundary::periodic)    return (f[l - config.offset * (config.ny - 1)] - f[l]) / config.hy;
			else if (b.type == MathBoundary::array_value) return (b.v[l % config.offset] - f[l]) / (0.5 * config.hy);
		}
		else {
			static_assert(false, "Invalid Side template parameter");
		}
		return 0; // дл€ подавлени€ предупреждений компил€тора
	}

	double dx2_border_forward(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l + 1]) / (pow(config.hx * 0.5, 2));
	};
	double dx2_border_back(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l - 1]) / (pow(config.hx * 0.5, 2));
	};
	double dy2_border_up(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l + config.offset]) / (pow(config.hy * 0.5, 2));
	};
	double dy2_border_down(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l - config.offset]) / (pow(config.hy * 0.5, 2));
	};

	double dx1_forward(double* f, unsigned int l, double borderF = 0)
	{
		return (f[l] - borderF) / (0.5 * config.hx);
	}
	double dx1_back(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - f[l]) / (0.5 * config.hx);
	}
	double dy1_up(double* f, unsigned int l, double borderF = 0)
	{
		return (f[l] - borderF) / (0.5 * config.hy);
	}
	double dy1_down(double* f, unsigned int l, double borderF = 0)
	{
		return (borderF - f[l]) / (0.5 * config.hy);
	}


	double dx1(Variable& f, unsigned int i, unsigned int l)
	{
		return (value2<Side::east>(f, l, i == config.nx - 1, f.wall[east]) - value2<Side::west>(f, l, i == 0, f.wall[west])) / (config.hx);
	}
	double dy1(Variable& f, unsigned int j, unsigned int l)
	{
		return (value2<Side::north>(f, l, j == config.ny - 1, f.wall[north]) - value2<Side::south>(f, l, j == 0, f.wall[south])) / (config.hy);
	}

	#define DF(side, wall) deriv<Side::side>(f, l, wall, f_bc[side].type, f_bc[side].value) 
	#define FF(side, wall) value<Side::side>(f, l, wall, f_bc[side].type, f_bc[side].value) 

	#define DER(f, side, if_boundary) deriv2<Side::side>(f, l, if_boundary, f.wall[side]) 
	#define VAL(f, side, if_boundary) value2<Side::side>(f, l, if_boundary, f.wall[side]) 

	#define DFl(side, wall, l) deriv<Side::side>(f, l, wall, f_bc[side].type, f_bc[side].value) 
	#define FFl(side, wall, l) value<Side::side>(f, l, wall, f_bc[side].type, f_bc[side].value) 

	#define DUX(side, wall) deriv<Side::side>(vx, l, wall, vx_bc[side].type, vx_bc[side].value) 
	#define FUX(side, wall) value<Side::side>(vx, l, wall, vx_bc[side].type, vx_bc[side].value)

	#define DUY(side, wall) deriv<Side::side>(vy, l, wall, vy_bc[side].type, vy_bc[side].value)
	#define FUY(side, wall) value<Side::side>(vy, l, wall, vy_bc[side].type, vy_bc[side].value)

	#define SX config.Sx
	#define SY config.Sy
	#define OFF config.offset


	void form_H_variable_viscosity(Configuration& c, Variable& H, Variable& vx, Variable& vy)
	{
		unsigned int l; // ll;
		for (unsigned int j = 0; j < c.ny; j++) {
			for (unsigned int i = 0; i < c.nx; i++) {
				l = i + c.offset * j;

				double tau0 = 250;
				double alpha = 5;
				double coef;

				coef = alpha + sqrt(2 * pow(dx1(vx, i, l), 2) + 2 * pow(dy1(vy, j, l), 2) + pow(dx1(vy, i, l) + dy1(vx, j, l), 2));
				H[l] = 1.0 + tau0 / coef;

				if (j == 0)
				{
					coef = alpha + sqrt(2 * pow(dy1_up(vy, l), 2) + pow(dy1_up(vx, l), 2));
					H.wall[south](i) = 1.0 + tau0 / coef;
				}
				if (j == c.ny - 1)
				{
					coef = alpha + sqrt(2 * pow(dy1_down(vy, l), 2) + pow(dy1_down(vx, l), 2));
					H.wall[north](i) = 1.0 + tau0 / coef;
				}
				if (i == 0)
				{
					coef = alpha + sqrt(2 * pow(dx1_forward(vx, l), 2) + pow(dx1_forward(vy, l), 2));
					H.wall[west](j) = 1.0 + tau0 / coef;
				}
				if (i == c.nx - 1)
				{
					coef = alpha + sqrt(2 * pow(dx1_back(vx, l), 2) + pow(dx1_back(vy, l), 2));
					H.wall[east](j) = 1.0 + tau0 / coef;
				}
			}
		}
	}


	void write_fields2d_with_boundaries(Configuration& c, std::string folder, std::string name, std::vector<Variable> v = {}, std::string head = "T, vx, vy, vz, C, Psi")
	{
		string str = folder + name + ".txt";

		std::ofstream all(str.c_str());

		all << "x, y, " + head << endl;
		//all << name << endl;
		//all << setprecision(16) << fixed;
		for (unsigned int j = 0; j < c.ny + c.q; j++) {
			for (unsigned int i = 0; i < c.nx + c.q; i++) {
				unsigned int l = i + c.offset * j;

				double x = x_coord(c, i);
				double y = y_coord(c, j);

				all << x << " " << y << " ";
				for (auto& it : v)
				{
					all << it[l] << " ";
				}
				all << endl;
			}
		}


		for (unsigned int j = 0; j < c.ny + c.q; j++)
		{
			double x = 0;
			double y = y_coord(c, j);
			all << x << " " << y << " ";
			for (auto& it : v)
				if (it.wall[west].type == MathBoundary::array_value)
					all << it.wall[west](j) << " ";
			all << endl;
		}

		for (unsigned int j = 0; j < c.ny + c.q; j++)
		{
			double x = c.Lx;
			double y = y_coord(c, j);
			all << x << " " << y << " ";

			for (auto& it : v)
				if (it.wall[east].type == MathBoundary::array_value)
					all << it.wall[east](j) << " ";
			all << endl;
		}

		for (unsigned int i = 0; i < c.nx + c.q; i++)
		{
			double x = x_coord(c, i);
			double y = 0;
			all << x << " " << y << " ";

			for (auto& it : v)
				if (it.wall[south].type == MathBoundary::array_value)
					all << it.wall[south](i) << " ";
			all << endl;
		}
		for (unsigned int i = 0; i < c.nx + c.q; i++)
		{
			double x = x_coord(c, i);
			double y = c.Ly;
			all << x << " " << y << " ";

			for (auto& it : v)
				if (it.wall[north].type == MathBoundary::array_value)
					all << it.wall[north](i) << " ";
			all << endl;
		}

	};



	double diffusion(unsigned int l, unsigned int i, unsigned int j, Variable& f)
	{
		double diff = 0.0;
		diff += SX * (DER(f, east, i == config.nx - 1) - DER(f, west, i == 0));
		diff += SY * (DER(f, north, j == config.ny - 1) - DER(f, south, j == 0));
		return diff;
	}

	double diffusion(unsigned int l, unsigned int i, unsigned int j, double* f, walls& f_bc)
	{
		double diff = 0.0;
		diff += (DF(east, i == config.nx - 1) - DF(west, i == 0)) * SX;
		diff += (DF(north, j == config.ny - 1) - DF(south, j == 0)) * SY;
		return diff;
	}
	double diffusion_variable_coef(unsigned int l, unsigned int i, unsigned int j, Variable& f, Variable& H)
	{
		double diff = 0.0;
		diff += SX * (deriv2<Side::east>(f, l, i == config.nx - 1, f.wall[east]) * value2<Side::east>(H, l, i == config.nx - 1, H.wall[east])
			- deriv2<Side::west>(f, l, i == 0, f.wall[west]) * value2<Side::west>(H, l, i == 0, H.wall[west]));
		diff += SY * (deriv2<Side::north>(f, l, j == config.ny - 1, f.wall[north]) * value2<Side::north>(H, l, j == config.ny - 1, H.wall[north])
			- deriv2<Side::south>(f, l, j == 0, f.wall[south]) * value2<Side::south>(H, l, j == 0, H.wall[south]));
		return diff;
	}
	double conv_central(unsigned int l, unsigned int i, unsigned int j, double* vx, double* vy, walls& vx_bc, walls& vy_bc, double* f, walls& f_bc)
	{
		double conv = 0.0;
		conv += (FF(east, i == config.nx - 1) * FUX(east, i == config.nx - 1) - FF(west, i == 0) * FUX(west, i == 0)) * SX;
		conv += (FF(north, j == config.ny - 1) * FUY(north, j == config.ny - 1) - FF(south, j == 0) * FUY(south, j == 0)) * SY;
		return conv;
	}
	double conv_central(unsigned int l, unsigned int i, unsigned int j, Variable& vx, Variable& vy, Variable& f)
	{
		double conv = 0.0;
		conv += SX * (value2<Side::east>(f, l, i == config.nx - 1, f.wall[east])
			* value2<Side::east>(vx, l, i == config.nx - 1, vx.wall[east])
			- value2<Side::west>(f, l, i == 0, f.wall[west])
			* value2<Side::west>(vx, l, i == 0, vx.wall[west]));
		conv += SY * (value2<Side::north>(f, l, j == config.ny - 1, f.wall[north])
			* value2<Side::north>(vy, l, j == config.ny - 1, vy.wall[north])
			- value2<Side::south>(f, l, j == 0, f.wall[south])
			* value2<Side::south>(vy, l, j == 0, vy.wall[south]));
		return conv;
	}

	double conv_upwind(unsigned int l, unsigned int i, unsigned int j, double* vx, double* vy, walls& vx_bc, walls& vy_bc, double* f, walls& f_bc)
	{
		double Fe, Fw, Fn, Fs;

		double Uxe = FUX(east, i == config.nx - 1);
		double Uxw = FUX(west, i == 0);

		double Uyn = FUY(north, j == config.ny - 1);
		double Uys = FUY(south, j == 0);

		if (Uxw > 0)
		{
			if (i != 0)	Fw = Uxw * f[l - 1];
			else		Fw = 0;
		}
		else
		{
			Fw = Uxw * f[l];
		}

		if (Uxe > 0)
		{
			if (i != config.nx - 1)	Fe = Uxe * f[l];
			else					Fe = 0;
		}
		else
		{
			Fe = Uxe * f[l + 1];
		}


		if (Uys > 0)
		{
			if (j != 0)	Fs = Uys * f[l - OFF];
			else		Fs = 0;
		}
		else
		{
			Fs = Uys * f[l];
		}

		if (Uyn > 0)
		{
			if (j != config.ny - 1)	Fn = Uyn * f[l];
			else					Fn = 0;
		}
		else
		{
			Fn = Uyn * f[l + OFF];
		}
		return (Fe - Fw) * SX + (Fn - Fs) * SY;
	}
	double div(unsigned int l, unsigned i, unsigned j, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
	{
		return (FUX(east, i == config.nx - 1) - FUX(west, i == 0)) * SX
			+ (FUY(north, j == config.ny - 1) - FUY(south, j == 0)) * SY;
	}
	void make_conv_buffer(Configuration& config, double* buffer, Variable& vx, Variable& vy, Variable& f)
	{
		for (unsigned int j = 0; j < config.ny; j++) {
			for (unsigned int i = 0; i < config.nx; i++) {
				unsigned int l = i + config.offset * j;
				buffer[l] = -conv_central(l, i, j, vx, vy, f) * config.tau / config.dV;
				//buffer[l] = -conv_upwind(l, i, j, vx, vy, vx_bc, vy_bc, f, f_bc) * config.tau / config.dV;
			}
		}
	}

	void swap_one(Configuration& config, double* fnew, double* fold)
	{
		for (unsigned int l = 0; l < config.N; l++)
		{
			fold[l] = fnew[l];
		}
	}
	void solve_jacobi(double* f, double* f0, double* bb, int NN, SparseMatrix& M)
	{
		unsigned int k = 0;
		const double eps_iter = 1e-5;

		for (k = 1; k < 10000; k++)
		{
			double s = 0;
			for (int j = 0; j < NN; j++)
			{
				s = M.line(j, f0);
				f[j] = f0[j] + (bb[j] - s) / M[j][j];
			}


			double max = 0;
			double dif;
			for (int i = 0; i < NN; i++)
			{
				dif = abs(f0[i] - f[i]);
				if (dif > max)
					max = dif;
			}
			for (int j = 0; j < NN; j++)
				f0[j] = f[j];

			if (max < eps_iter)	break;

			if (k % 1000 == 0) cout << "host k = " << k << ", eps = " << max << endl;
		}

	}


	namespace projection_convection
	{
		void temperature(Configuration& config, Variable& T_new, Variable& T, Variable& vx, Variable& vy)
		{
			double coef = 1.0 / config.Pr;
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;
					double diff = diffusion(l, i, j, T);
					double conv = conv_central(l, i, j, vx, vy, T);

					T_new[l] = T[l] + config.tau / config.dV * (diff / config.Pr - conv);
				}
			}
		}
		void quasi_velocity(Configuration& config, Variable& ux, Variable& uy, Variable& vx, Variable& vy, Variable& T)
		{
			double diff = 0, conv = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					diff = diffusion(l, i, j, vx);
					conv = conv_central(l, i, j, vx, vy, vx);
					ux[l] = vx[l] + config.tau / config.dV * (diff - conv);

					diff = diffusion(l, i, j, vy);
					conv = conv_central(l, i, j, vx, vy, vy);
					uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.Ra / config.Pr * config.tau * T[l];
				}
			}
		}
		void quasi_velocity_variable_viscosity(Configuration& config, Variable& ux, Variable& uy, Variable& vx, Variable& vy, Variable& T, Variable& H)
		{
			double diff = 0, conv = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					diff = diffusion_variable_coef(l, i, j, vx, H);
					//diff = diffusion(l, i, j, vx);
					conv = conv_central(l, i, j, vx, vy, vx);
					ux[l] = vx[l] + config.tau / config.dV * (diff - conv);

					diff = diffusion_variable_coef(l, i, j, vy, H);
					//diff = diffusion(l, i, j, vy);
					conv = conv_central(l, i, j, vx, vy, vy);
					uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.Ra / config.Pr * config.tau * T[l];
				}
			}
		}
		void poisson_pressure(Configuration& config, Variable& f_new, Variable& f, Variable& ux, Variable& uy)
		{
			config.tau_p = 0.25 * config.hx * config.hx;
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					double d = 0.0;

					if (i == 0)
					{
						if (config.domain == closed_box) d += 0.0;
						if (config.domain == open_tube)  d += -SX * (deriv2<Side::west>(f, l, true, f.wall[west]) - ux[l] / config.tau);
					}
					else
						d += -(deriv2<Side::west>(f, l) - value2<Side::west>(ux, l) / config.tau) * config.Sx;

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box) d += 0.0;
						if (config.domain == open_tube)  d += SX * (deriv2<Side::east>(f, l, true, f.wall[east]) - ux[l] / config.tau);
					}
					else
						d += +(deriv2<Side::east>(f, l) - value2<Side::east>(ux, l) / config.tau) * config.Sx;

					if (j == 0)						d += 0.0;
					else							d += -SY * (deriv2<Side::south>(f, l) - value2<Side::south>(uy, l) / config.tau);

					if (j == config.ny - 1)			d += 0.0;
					else							d += +SY * (deriv2<Side::north>(f, l) - value2<Side::north>(uy, l) / config.tau);


					//d = dPflux(f, l, i, j, ux, uy);

					f_new[l] = f[l] + config.tau_p * ((d) / config.dV);
				}
			}
		}
		void solve_poisson(Configuration& config, Variable& ksi, Variable& ksi0, Variable& ux, Variable& uy)
		{
			unsigned int k = 0;
			double eps = 0, res = 0, res0 = 0;
			double eps_iter = 1e-5;

			k = 0;
			eps = 1.0;
			res = 0.0;
			res0 = 0.0;

			auto reduce = [&](double* f)
			{
				double s = 0;
				for (unsigned int l = 0; l < config.N; l++)
					s += abs(f[l]);
				return s;
			};


			for (k = 1; k < 100000; k++)
				//while(true)
			{
				poisson_pressure(config, ksi, ksi, ux, uy);
				res = reduce(ksi);
				eps = abs(res - res0) / (res0 + 1e-6);
				res0 = res;

				//std::swap(ksi, ksi0);
				for (unsigned int l = 0; l < config.N; l++)
					ksi0[l] = ksi[l];

				if (eps < eps_iter)	break;
				if (k % 1000 == 0) std::cout << "Poisson, k = " << k << ", eps = " << eps << std::endl;
			}
			//if (k > 100) std::cout << "device k = " << k << ", eps = " << eps << std::endl;
		}
		void velocity_correction(Configuration& config, Variable& p, Variable& vx, Variable& vy, Variable& ux, Variable& uy, Variable& T)
		{
			auto dx2_border_forward = [&config](double* f, double borderF, unsigned int l)
			{
				return (borderF - 1.5 * f[l] + 0.5 * f[l + 1]) / (pow(config.hx * 0.5, 2)); //h = half step, should be
			};
			auto dx2_border_back = [&config](double* f, double borderF, unsigned int l)
			{
				return (borderF - 1.5 * f[l] + 0.5 * f[l - 1]) / (pow(config.hx * 0.5, 2)); //h = half step, should be
			};
			auto dy2_border_up = [&config](double* f, double borderF, unsigned int l)
			{
				return (borderF - 1.5 * f[l] + 0.5 * f[l + config.offset]) / (pow(config.hy * 0.5, 2)); //h = half step, should be
			};
			auto dy2_border_down = [&config](double* f, double borderF, unsigned int l)
			{
				return (borderF - 1.5 * f[l] + 0.5 * f[l - config.offset]) / (pow(config.hy * 0.5, 2)); //h = half step, should be
			};


			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l] - 0.5 * config.hx * dx2_border_forward(vx, 0, l);
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l] + config.hx * 0.5 * dx2_border_back(vx, 0, l);
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value<Side::east>(p, l);
					}


					if (j == 0)
					{
						//if (config.domain == closed_box) 
						Ps = p[l] - 0.5 * config.hy * (dy2_border_up(vy, 0, l) + config.Ra / config.Pr * T.wall[south].value);
					}
					else
					{
						Ps = value<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						//if (config.domain == closed_box) 
						Pn = p[l] + 0.5 * config.hy * (dy2_border_down(vy, 0, l) + config.Ra / config.Pr * T.wall[north].value);
					}
					else
					{
						Pn = value<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / config.dV;
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / config.dV;
				}
			}
		}
		void velocity_correction2(Configuration& config, Variable& p, Variable& vx, Variable& vy, Variable& ux, Variable& uy, Variable& T)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l] - 0.5 * config.hx * dx2_border_forward(vx, l);
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value2<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l] + config.hx * 0.5 * dx2_border_back(vx, l);
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value2<Side::east>(p, l);
					}

					if (j == 0)
					{
						//if (config.domain == closed_box) 
						Ps = p[l] - 0.5 * config.hy * (dy2_border_up(vy, l) + config.Ra / config.Pr * T.wall[south](i));
					}
					else
					{
						Ps = value2<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						//if (config.domain == closed_box) 
						Pn = p[l] + 0.5 * config.hy * (dy2_border_down(vy, l) + config.Ra / config.Pr * T.wall[north](i));
					}
					else
					{
						Pn = value2<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / config.dV;
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / config.dV;
				}
			}
		}

		void velocity_correction3(Configuration& config, Variable& p, Variable& vx, Variable& vy, Variable& ux, Variable& uy, Variable& T)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l]; // -0.5 * config.hx * dx2_border_forward(vx, l);
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value2<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l]; // +config.hx * 0.5 * dx2_border_back(vx, l);
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value2<Side::east>(p, l);
					}

					if (j == 0)
					{
						//if (config.domain == closed_box) 
						Ps = p[l]; // -0.5 * config.hy * (dy2_border_up(vy, l) + config.Ra / config.Pr * T.wall[south](i));
					}
					else
					{
						Ps = value2<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						//if (config.domain == closed_box) 
						Pn = p[l]; // +0.5 * config.hy * (dy2_border_down(vy, l) + config.Ra / config.Pr * T.wall[north](i));
					}
					else
					{
						Pn = value2<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / (config.dV * 0.5);
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / (config.dV * 0.5);
				}
			}
		}
		void quasi_velocity_fv_v2(Configuration& config, double* ux, double* uy, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
		{
			double diff = 0, conv = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					//diff = diffusion_variable_coef(l, i, j, vx, vx_bc) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vx_bc, vy_bc, vx, vx_bc);
					ux[l] = (diff - conv) / config.dV;
					ux[l] = ux[l] * config.tau + vx[l];

					diff = diffusion(l, i, j, vy, vy_bc) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vx_bc, vy_bc, vy, vy_bc);
					uy[l] = (diff - conv) / config.dV;
					uy[l] = uy[l] * config.tau + vy[l];
				}
			}
		}
	}

	namespace projection_phase
	{
		inline double phase_src_x_(Configuration& c, unsigned int i, unsigned int l, double* C, double* mu)
		{
			double mu1 = (i == 0) ? mu[l] : mu[l - 1];
			double mu2 = (i == config.nx - 1) ? mu[l] : mu[l + 1];

			return -C[l] * (mu2 - mu1) / (config.hx * 2) / config.Mach;
		}
		inline double phase_src_y_(Configuration& c, unsigned int j, unsigned int l, double* C, double* mu)
		{
			double mu1 = (j == 0) ? mu[l] : mu[l - c.offset];
			double mu2 = (j == config.ny - 1) ? mu[l] : mu[l + c.offset];

			return -C[l] * (mu2 - mu1) / (config.hy * 2) / config.Mach + config.Gr * C[l];
		}
		inline double phase_src_x(Configuration& c, unsigned int i, unsigned int l, Variable& C, Variable& mu)
		{
			return -C[l] * (VAL(mu, east, i == c.nx - 1) - VAL(mu, west, i == 0)) / (c.hx) / c.Mach;
		}
		inline double phase_src_y(Configuration& c, unsigned int j, unsigned int l, Variable& C, Variable& mu)
		{
			return -C[l] * (VAL(mu, north, j == c.ny - 1) - VAL(mu, south, j == 0)) / (c.hy) / c.Mach + c.Gr * C[l];
		}

		inline double phase_src_x_cut(Configuration& c, unsigned int i, unsigned int l, Variable& C, Variable& mu)
		{
			double coef = 1;
			if (i == c.nx - 1 || i == 0) coef = 0.5;
			return -C[l] * (VAL(mu, east, i == c.nx - 1) - VAL(mu, west, i == 0)) / (c.hx * coef) / c.Mach;
		}
		inline double phase_src_y_cut(Configuration& c, unsigned int j, unsigned int l, Variable& C, Variable& mu)
		{
			double gr = c.Gr * C[l];
			double coef = 1;
			if (j == c.ny - 1 || j == 0) coef = 0.5;
			return -C[l] * (VAL(mu, north, j == c.ny - 1) - VAL(mu, south, j == 0)) / (c.hy * coef) / c.Mach + gr;
		}


		void form_H_viscosity(Configuration& c, double* H, double* C, double* buffer = nullptr, double* buffer2 = nullptr)
		{
			auto H_tanh_mirrored = [](double x, double Amin = 0.0, double Amax = 1.0, double x_shift = 0.5, double eps = 0.01)
			{
				double tanh_val = tanh((x - x_shift) / eps);
				return Amin + (Amax - Amin) * (1.0 - tanh_val) / 2.0;
			};

			unsigned int l; // ll;
			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					l = i + c.offset * j;
					H[l] = H_tanh_mirrored(C[l], 1, 1, 0.5, 0.1);
				}
			}
		}

		void make_source_for_phase_field(Configuration& config, double* src_x, double* src_y, double* C, double* mu)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;
					src_x[l] = phase_src_x_(config, i, l, C, mu);
					src_y[l] = phase_src_y_(config, j, l, C, mu);
				}
			}
		}
		void quasi_velocity(Configuration& config, Variable& ux, Variable& uy, Variable& vx, Variable& vy, Variable& C, Variable& mu)
		{
			double diff = 0, conv = 0, src_x = 0, src_y = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					diff = diffusion(l, i, j, vx) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vx);
					src_x = phase_src_x_cut(config, i, l, C, mu);
					//src_x = -C[l] * (VAL(mu, east, i == config.nx - 1) - VAL(mu, west, i == 0)) / (config.hx) / config.Mach;
					ux[l] = vx[l] + config.tau / config.dV * (diff - conv) + config.tau * src_x;

					diff = diffusion(l, i, j, vy) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vy);
					src_y = phase_src_y_cut(config, j, l, C, mu);
					//src_y = -C[l] * (VAL(mu, north, j == config.ny - 1) - VAL(mu, south, j == 0)) / (config.hy) / config.Mach;
					//src_y += config.Gr * C[l];
					uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.tau * src_y;
				}
			}
		}

		inline double phase3_src_x_cut(Configuration& c, unsigned int i, unsigned int l, Variable& C, Variable& mu, Variable& C2, Variable& mu2, Variable& C3, Variable& mu3)
		{
			double coef = 1;
			if (i == c.nx - 1 || i == 0) coef = 0.5;
			double val = 0;
			val += -C[l] * (VAL(mu, east, i == c.nx - 1) - VAL(mu, west, i == 0)) / (c.hx * coef) / c.Mach;
			val += -C2[l] * (VAL(mu2, east, i == c.nx - 1) - VAL(mu2, west, i == 0)) / (c.hx * coef) / c.Mach2;
			val += -C3[l] * (VAL(mu3, east, i == c.nx - 1) - VAL(mu3, west, i == 0)) / (c.hx * coef) / c.Mach3;
			return val;
		}
		inline double phase3_src_y_cut(Configuration& c, unsigned int j, unsigned int l, Variable& C, Variable& mu, Variable& C2, Variable& mu2, Variable& C3, Variable& mu3)
		{
			double gr = c.Gr * C[l] + c.Gr2 * C2[l] + c.Gr3 * C3[l];
			double coef = 1;
			if (j == c.ny - 1 || j == 0) coef = 0.5;
			double val = gr;
			val += -C[l] * (VAL(mu, north, j == c.ny - 1) - VAL(mu, south, j == 0)) / (c.hy * coef) / c.Mach;
			val += -C2[l] * (VAL(mu2, north, j == c.ny - 1) - VAL(mu2, south, j == 0)) / (c.hy * coef) / c.Mach2;
			val += -C3[l] * (VAL(mu3, north, j == c.ny - 1) - VAL(mu3, south, j == 0)) / (c.hy * coef) / c.Mach3;
			return val;
		}
		void quasi_velocity_3_phases(Configuration& config, Variable& ux, Variable& uy, Variable& vx, Variable& vy, Variable& C, Variable& mu, Variable& C2, Variable& mu2, Variable& C3, Variable& mu3)
		{
			double diff = 0, conv = 0, src_x = 0, src_y = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					diff = diffusion(l, i, j, vx) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vx);
					src_x = phase3_src_x_cut(config, i, l, C, mu, C2, mu2, C3, mu3);
					//src_x = -C[l] * (VAL(mu, east, i == config.nx - 1) - VAL(mu, west, i == 0)) / (config.hx) / config.Mach;
					ux[l] = vx[l] + config.tau / config.dV * (diff - conv) + config.tau * src_x;

					diff = diffusion(l, i, j, vy) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vy);
					src_y = phase3_src_y_cut(config, j, l, C, mu, C2, mu, C3, mu3);
					//src_y = -C[l] * (VAL(mu, north, j == config.ny - 1) - VAL(mu, south, j == 0)) / (config.hy) / config.Mach;
					//src_y += config.Gr * C[l];
					uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.tau * src_y;
				}
			}
		}


		void quasi_velocity_H(Configuration& config, Variable& ux, Variable& uy, Variable& vx, Variable& vy, Variable& C, Variable& mu, Variable& H)
		{
			double diff = 0, conv = 0, src_x = 0, src_y = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					//diff = diffusion(l, i, j, vx) / config.Re;
					diff = diffusion_variable_coef(l, i, j, vx, H) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vx);
					src_x = phase_src_x(config, i, l, C, mu);
					ux[l] = vx[l] + config.tau / config.dV * (diff - conv) + config.tau * src_x;

					//diff = diffusion(l, i, j, vy) / config.Re;
					diff = diffusion_variable_coef(l, i, j, vy, H) / config.Re;
					conv = conv_central(l, i, j, vx, vy, vy);
					src_y = phase_src_y(config, j, l, C, mu);
					uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.tau * src_y;
				}
			}
		}


		void velocity_correction_phase(Configuration& config, double* p, double* vx, double* vy, double* ux, double* uy, double* f, walls& C_bc)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l] - 0.5 * config.hx * dx2_border_forward(vx, l, 0);
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l] + config.hx * 0.5 * dx2_border_back(vx, l, 0);
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value<Side::east>(p, l);
					}


					if (j == 0)
					{
						//if (config.domain == closed_box) 
						Ps = p[l] - 0.5 * config.hy * dy2_border_up(vy, l, 0);
						Ps += -config.Gr * 0.5 * config.hy * value<Side::south>(f, l, true, C_bc[south].type, C_bc[south].value);
					}
					else
					{
						Ps = value<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						//if (config.domain == closed_box) 
						Pn = p[l] + 0.5 * config.hy * dy2_border_down(vy, l, 0);
						Pn += config.Gr * 0.5 * config.hy * value<Side::north>(f, l, true, C_bc[north].type, C_bc[north].value);
					}
					else
					{
						Pn = value<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / config.dV;
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / config.dV;
				}
			}
		}
		void velocity_correction_cut(Configuration& config, double* p, double* vx, double* vy, double* ux, double* uy, Variable& C)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l];
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l];
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value<Side::east>(p, l);
					}


					if (j == 0)
					{
						Ps = p[l];
					}
					else
					{
						Ps = value<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						Pn = p[l];
					}
					else
					{
						Pn = value<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / (config.dV * 0.5);
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / (config.dV * 0.5);
				}
			}
		}

		void velocity_correction_phase3(Configuration& config, double* p, double* vx, double* vy, double* ux, double* uy, Variable& C, Variable& C2)
		{
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					unsigned int l = i + config.offset * j;
					double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

					if (i == 0)
					{
						if (config.domain == closed_box) Pw = p[l] - 0.5 * config.hx * dx2_border_forward(vx, l, 0);
						if (config.domain == open_tube)  Pw = config.p_in;
					}
					else
					{
						Pw = value<Side::west>(p, l);
					}

					if (i == config.nx - 1)
					{
						if (config.domain == closed_box)  Pe = p[l] + config.hx * 0.5 * dx2_border_back(vx, l, 0);
						if (config.domain == open_tube)	  Pe = config.p_out;
					}
					else
					{
						Pe = value<Side::east>(p, l);
					}


					if (j == 0)
					{
						//if (config.domain == closed_box) 
						Ps = p[l] - 0.5 * config.hy * dy2_border_up(vy, l, 0);
						Ps += -config.Gr * 0.5 * config.hy * value2<Side::south>(C, l, true, C.wall[south]);
						Ps += -config.Gr2 * 0.5 * config.hy * value2<Side::south>(C2, l, true, C2.wall[south]);
					}
					else
					{
						Ps = value<Side::south>(p, l);
					}

					if (j == config.ny - 1)
					{
						//if (config.domain == closed_box) 
						Pn = p[l] + 0.5 * config.hy * dy2_border_down(vy, l, 0);
						Pn += config.Gr * 0.5 * config.hy * value2<Side::north>(C, l, true, C.wall[north]);
						Pn += config.Gr2 * 0.5 * config.hy * value2<Side::north>(C2, l, true, C2.wall[north]);
					}
					else
					{
						Pn = value<Side::north>(p, l);
					}

					vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / config.dV;
					vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / config.dV;
				}
			}
		}
	}

	namespace flow_test
	{
		void quasi_velocity_test(Configuration& config, double* ux_buf, double* uy_buf, double* ux, double* uy, double* vx, double* vy, double* f, walls& f_bc, walls& vx_bc, walls& vy_bc)
		{
			double diff = 0, conv = 0;

			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					diff = diffusion(l, i, j, ux, vx_bc);
					conv = conv_upwind(l, i, j, ux, uy, vx_bc, vy_bc, ux, vx_bc);
					ux_buf[l] = vx[l] + config.tau / config.dV * (diff - conv);

					diff = diffusion(l, i, j, uy, vy_bc);
					conv = conv_upwind(l, i, j, ux, uy, vx_bc, vy_bc, uy, vy_bc);
					uy_buf[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.Ra / config.Pr * config.tau * f[l];
				}
			}
		}
		void form_rhs_for_heat_equation(Configuration& c, double* b, double ax, double ay, double* vx, double* vy, walls& vx_bc, walls& vy_bc, double* main, walls& main_bc)
		{
			unsigned int l;
			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					l = i + c.offset * j;
					b[l] = 0;

					auto bc = [&](enumside side, double a, bool border, walls& w)
					{
						if (border)
						{
							if (w[side].type == MathBoundary::Neumann)
							{
								b[l] += w[side].value;
							}
							else if (w[side].type == MathBoundary::Dirichlet)
							{
								b[l] += w[side].value * (a / 0.5);
							}
							else if (w[side].type == MathBoundary::periodic)
							{

							}
						}
					};

					bc(west, ax, i == 0, main_bc);
					bc(east, ax, i == c.nx - 1, main_bc);

					bc(south, ay, j == 0, main_bc);
					bc(north, ay, j == c.ny - 1, main_bc);

					b[l] += c.dV / c.tau * main[l];
					b[l] += -conv_central(l, i, j, vx, vy, vx_bc, vy_bc, main, main_bc);
				}
			}
		}

		//  Scheme:	(dV / tau + Ap) * U;  Ae = Aw = - S / h * coef;  Ap = + 2 * S / h * coef 
		void form_matrix_with_2d_laplace_FV(Configuration& c, SparseMatrix& M, /*const char* str,*/ double ax, double ay, double at, bc_type bcx, bc_type bcy)
		{
			unsigned int l;
			//int type = 0;
			//if (str == "diagonal") type = 0;
			//if (str == "unity")	   type = 1;
			auto index = [&](Side side, unsigned int l, bool periodic = false)
			{
				if (side == Side::center)
					return l;

				else if (!periodic)
				{
					switch (side)
					{
						//case Side::center:
						//	return l;
						//	break;
					case Side::west:
						return l - 1;
						break;
					case Side::east:
						return l + 1;
						break;
					case Side::south:
						return l - c.offset;
						break;
					case Side::north:
						return l + c.offset;
						break;
					case Side::front:
						return l - c.offset2;
						break;
					case Side::back:
						return l + c.offset2;
						break;
					}
				}
				else if (periodic)
				{
					switch (side)
					{
						//case Side::center:
						//	return l;
						//	break;
					case Side::west:
						return l + (c.nx - 1);
						break;
					case Side::east:
						return l - (c.nx - 1);
						break;
					case Side::south:
						return l + (c.ny - 1) * c.offset;
						break;
					case Side::north:
						return l - (c.ny - 1) * c.offset;
						break;
					case Side::front:
						return l + (c.nz - 1) * c.offset2;
						break;
					case Side::back:
						return l - (c.nz - 1) * c.offset2;
						break;
					}
				}
				return UINT_MAX;
			};

			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					l = i + c.offset * j;
					std::map<int, double> m;
					auto flux = [&](Side side, double a, bool border, bc_type bc)
					{
						if (border)
						{
							if (bc == Neumann)
							{
								//
							}
							if (bc == Dirichlet)
							{
								m[index(Side::center, l)] += a / 0.5;
							}
							if (bc == periodic)
							{
								m[index(Side::center, l)] += a;
								m[index(side, l, true)] = -a;
							}
						}
						else
						{
							m[index(Side::center, l)] += a;
							m[index(side, l)] = -a;
						}
					};

					flux(Side::west, ax, i == 0, bcx);
					flux(Side::east, ax, i == c.nx - 1, bcx);

					flux(Side::south, ay, j == 0, bcy);
					flux(Side::north, ay, j == c.ny - 1, bcy);

					m[index(Side::center, l)] += at;


					M.add_line_with_map(m, l);
				}
				if (j % 10 == 0 || j == c.ny - 1)
					cout << "matrix ready: " << 100.0 * (j + 1) / c.ny << "%" << endl;
			}
		}
		void form_matrix_with_2d_laplace_FV(Configuration& c, SparseMatrix& M, /*const char* str,*/ double ax, double ay, double at, walls& bc)
		{
			unsigned int l;
			//int type = 0;
			//if (str == "diagonal") type = 0;
			//if (str == "unity")	   type = 1;

			auto index = [&](enumside side, unsigned int l, bool periodic = false)
			{
				if (side == inner)
					return l;

				else if (!periodic)
				{
					switch (side)
					{
						//case Side::center:
						//	return l;
						//	break;
					case west:
						return l - 1;
						break;
					case east:
						return l + 1;
						break;
					case south:
						return l - c.offset;
						break;
					case north:
						return l + c.offset;
						break;
					case front:
						return l - c.offset2;
						break;
					case back:
						return l + c.offset2;
						break;
					}
				}
				else if (periodic)
				{
					switch (side)
					{
						//case Side::center:
						//	return l;
						//	break;
					case west:
						return l + (c.nx - 1);
						break;
					case east:
						return l - (c.nx - 1);
						break;
					case south:
						return l + (c.ny - 1) * c.offset;
						break;
					case north:
						return l - (c.ny - 1) * c.offset;
						break;
					case front:
						return l + (c.nz - 1) * c.offset2;
						break;
					case back:
						return l - (c.nz - 1) * c.offset2;
						break;
					}
				}
				return UINT_MAX;
			};

			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					l = i + c.offset * j;
					std::map<int, double> m;


					auto flux = [&](enumside side, double a, bool border, walls& w)
					{
						if (border)
						{
							if (w[side].type == MathBoundary::Neumann)
							{
								//
							}
							if (w[side].type == MathBoundary::Dirichlet)
							{
								m[index(inner, l)] += a / 0.5;
							}
							if (w[side].type == MathBoundary::periodic)
							{
								m[index(inner, l)] += a;
								m[index(side, l, true)] = -a;
							}
						}
						else
						{
							m[index(inner, l)] += a;
							m[index(side, l)] = -a;
						}
					};

					flux(west, ax, i == 0, bc);
					flux(east, ax, i == c.nx - 1, bc);

					flux(south, ay, j == 0, bc);
					flux(north, ay, j == c.ny - 1, bc);

					m[index(inner, l)] += at;
					M.add_line_with_map(m, l);
				}
				if (j % 10 == 0 || j == c.ny - 1)
					cout << "matrix ready: " << 100.0 * (j + 1) / c.ny << "%" << endl;
			}
		}

		void filter(Configuration& c, double* p_new, double* p, walls& p_bc)
		{
			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++)
				{
					unsigned int l = i + c.offset * j;
					//double pe = DER(p, east, i == c.nx - 1);
					//double pw = DER(p, west, i == 0);
					//double pn = DER(p, north, j == c.ny - 1);
					//double ps = DER(p, south, j == 0);

					//p_new[l] = p[l] + ((pe - pw) * SX + (pn - ps) * SY) / c.dV * 0.001;
				}
			}
		}
		double check_div(Configuration& c, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
		{
			double max = 0;
			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					unsigned int l = i + c.offset * j;
					double d = div(l, i, j, vx, vy, vx_bc, vy_bc);
					if (abs(d) > max) max = abs(d);
				}
			}
			return max;
		}

	}

	struct FlowSolver : Configuration, Arrays
	{
		walls T_bc, vx_bc, vy_bc, p_bc, C_bc, mu_bc, p_prime_bc, vx_prime_bc, vy_prime_bc;
		SparseMatrix mFlow;
		SparseMatrix mTest;
		SparseMatrix pSM;
		bool ready = false;

		FlowSolver() : Configuration() {}
		FlowSolver(const Configuration& config, const Arrays& arr) : Configuration(config), Arrays(arr)
		{
			init_bc({ &T_bc, &vx_bc, &vy_bc, &p_bc, &C_bc, &mu_bc, &p_prime_bc, &vx_prime_bc, &vy_prime_bc });

			pSM.resize(N);
			double alpha = 0.02;
			flow_test::form_matrix_with_2d_laplace_FV(*this, pSM, alpha, alpha, 1, p_prime_bc);
		}
		void make_flow_matrix(bool load = false)
		{
			SparseMatrix mUX(N), mUY(N);

			double ax = (Sx / hx) / Re;
			double ay = (Sy / hy) / Re;
			double at = dV / tau;

			if (load == true)
			{
				mFlow.read_compressed_matrix("SparseMatrix_NS.dat");
			}

			if (load == false)
			{
				if (domain == closed_box)
				{
					//form_matrix_with_2d_laplace_FV(*this, mUX, ax, ay, at, Dirichlet, Dirichlet);
					//form_matrix_with_2d_laplace_FV(*this, mUY, ax, ay, at, Dirichlet, Dirichlet);
					flow_test::form_matrix_with_2d_laplace_FV(*this, mUX, ax, ay, at, vx_bc);
					flow_test::form_matrix_with_2d_laplace_FV(*this, mUY, ax, ay, at, vy_bc);

					mFlow.make_sparse_from_joint({ &mUX, &mUY });
				}
				else
				{
					cout << "Error: not yet implemented " << endl;
				}

				mFlow.save_compressed_matrix("SparseMatrix_NS.dat");
			}
		}
		void form_rhs_flow(Configuration& c, Component comp, unsigned int stride, double* b, double ax, double ay, double* vx, double* vy, walls& vx_bc, walls& vy_bc, double* p, walls& p_bc, double* vmain, walls& vmain_bc, double* T, walls& T_bc)
		{
			unsigned int l;
			for (unsigned int j = 0; j < c.ny; j++) {
				for (unsigned int i = 0; i < c.nx; i++) {
					l = i + c.offset * j + stride;
					b[l] = 0;

					auto bc = [&](enumside side, double a, bool border, walls& w)
					{
						if (border)
						{
							if (w[side].type == MathBoundary::Neumann)
							{
								b[l] += w[side].value;
							}
							else if (w[side].type == MathBoundary::Dirichlet)
							{
								b[l] += 0;// w[side].value* (a / 0.5);
							}
							else if (w[side].type == MathBoundary::periodic)
							{

							}
						}
					};

					bc(west, ax, i == 0, vmain_bc);
					bc(east, ax, i == c.nx - 1, vmain_bc);

					bc(south, ay, j == 0, vmain_bc);
					bc(north, ay, j == c.ny - 1, vmain_bc);

					b[l] += c.dV / c.tau * vmain[l];
					b[l] += -conv_central(l, i, j, vx, vy, vx_bc, vy_bc, vmain, vmain_bc);
					//b[l] += -conv_upwind(l, i, j, vx, vy, vx_bc, vy_bc, vmain, vmain_bc);

					if (comp == Component::x)
					{
						//b[l] += -(VAL(p, east, i == c.nx - 1) - VAL(p, west, i == 0)) * SX;
						//b[l] += phase_src_x(config, i, l, C, mu);
					}
					else if (comp == Component::y)
					{
						//b[l] += -(VAL(p, north, j == c.ny - 1) - VAL(p, south, j == 0)) * SY + c.Ra / c.Pr * T[l];
						//b[l] += phase_src_y(config, j, l, C, mu);
					}
				}
			}
		}
		void poisson_for_p_prime(Configuration& config, double* f_new, double* f, walls& f_bc, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
		{
			config.tau_p = 0.2 * config.hx * config.hx;
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++) {
					unsigned int l = i + config.offset * j;

					double d = 0.0;

					//if (i == 0)
					//{
					//	if (config.domain == closed_box) d += f[l];
					//	if (config.domain == open_tube)  d += 0.0;
					//}
					//else
					//	d += -(deriv<Side::west>(f, l) - value<Side::west>(ux, l) / config.tau) * config.Sx;

					//if (i == config.nx - 1)
					//{
					//	if (config.domain == closed_box) d += 0.0;
					//	if (config.domain == open_tube)  d += SX * (deriv<Side::east>(f, l, true, MathBoundary::Dirichlet, config.p_out) - ux[l] / config.tau);
					//}
					//else
					//	d += +(deriv<Side::east>(f, l) - value<Side::east>(ux, l) / config.tau) * config.Sx;

					//if (j == 0)						d += 0.0;
					//else							d += -SY * (deriv<Side::south>(f, l) - value<Side::south>(uy, l) / config.tau);

					//if (j == config.ny - 1)			d += 0.0;
					//else							d += +SY * (deriv<Side::north>(f, l) - value<Side::north>(uy, l) / config.tau);


					//f_new[l] = f[l] + config.tau_p * ((d) / config.dV);
					d = diffusion(l, i, j, f, f_bc) - div(l, i, j, vx, vy, vx_bc, vy_bc) / config.tau;
					f_new[l] = f[l] + config.tau_p / config.dV * (d);
				}
			}
		}
		void solve_poisson_for_p_prime(Configuration& config, double* p_new, double* p, walls& p_bc, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
		{
			unsigned int k = 0;
			double eps = 0, res = 0, res0 = 0;
			const double eps_iter = 1e-5;

			k = 0;
			eps = 1.0;
			res = 0.0;
			res0 = 0.0;

			auto reduce = [&](double* f)
			{
				double s = 0;
				for (unsigned int l = 0; l < config.N; l++)
					s += abs(f[l]);
				return s;
			};


			for (k = 1; k < 100000; k++)
				//while(true)
			{
				poisson_for_p_prime(config, p_new, p_new, p_bc, vx, vy, vx_bc, vy_bc);
				res = reduce(p_new);
				eps = abs(res - res0) / (res0 + 1e-5);
				res0 = res;

				//std::swap(ksi, ksi0);
				for (unsigned int l = 0; l < config.N; l++)
					p[l] = p_new[l];

				if (eps < eps_iter)	break;
				//if (k % 1000 == 0) std::cout << "k = " << k << ", eps = " << eps << std::endl;


			}
			//if (k > 100) std::cout << "device k = " << k << ", eps = " << eps << std::endl;
		}
		void correction(Configuration& config, double* p, walls& p_bc, double* p_prime, walls& p_prime_bc, double* vx, double* vy, walls& vx_bc, walls& vy_bc)
		{
			unsigned int l;
			double alpha = 0.05;
			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					l = i + config.offset * j;
					p[l] = p[l] + alpha * p_prime[l];
				}
			}

			//pressure_smoothing(p, p, p_prime);


			for (unsigned int j = 0; j < config.ny; j++) {
				for (unsigned int i = 0; i < config.nx; i++)
				{
					l = i + config.offset * j;

					//double pe = VAL(p_prime, east, i == config.nx - 1);
					//double pw = VAL(p_prime, west, i == 0);

					//double pn = VAL(p_prime, north, j == config.ny - 1);
					//double ps = VAL(p_prime, south, j == 0);

					//vx[l] = vx[l] - config.tau * (pe - pw) / config.hx;
					//vy[l] = vy[l] - config.tau * (pn - ps) / config.hy;
				}
			}
		}

		void pressure_smoothing(double* p_new, double* p0, double* p_prime)
		{
			for (unsigned int l = 0; l < N; l++)
				buffer[l] = p0[l] + p_prime[l];
			solve_jacobi(p_new, p_new, buffer, N, pSM);
		}
		void solve_simple()
		{
			unsigned int steps = 0;
			while (true)
			{
				steps++;

				//quasi_velocity_test(*this, ux_buf, uy_buf, ux, uy, vx, vy, T, T_bc, vx_bc, vy_bc);
				//quasi_velocity(*this, ux_buf, uy_buf, vx, vy, T, T_bc, vx_bc, vy_bc);
			   // solve_poisson(*this, p, p0, ux_buf, uy_buf);


				form_rhs_flow(*this, Component::x, 0, &Urhs[0], 1.0 / Re, 1.0 / Re, vx, vy, vx_bc, vy_bc, p, p_bc, vx, vx_bc, T, T_bc);
				form_rhs_flow(*this, Component::y, 0, &Urhs[N], 1.0 / Re, 1.0 / Re, vx, vy, vx_bc, vy_bc, p, p_bc, vy, vy_bc, T, T_bc);
				solve_jacobi(U, U, Urhs, N * 2, mFlow);

				for (unsigned int l = 0; l < N; l++) p_prime[l] = p_prime0[l] = 0.0;

				solve_poisson_for_p_prime(*this, p_prime, p_prime0, p_prime_bc, ux, uy, vx_bc, vy_bc);


				correction(*this, p, p_prime_bc, p_prime, p_prime_bc, ux, uy, vx_bc, vy_bc);

				double div_ = flow_test::check_div(*this, ux, uy, vx_bc, vy_bc);



				if (steps % 1 == 0)
				{
					//pause
					cout << "step = " << steps << ", div = " << div_ << endl;

					for (unsigned int j = 0; j < ny; j++) {
						for (unsigned int i = 0; i < nx; i++)
						{
							unsigned int l = i + offset * j;
							buffer[l] = div(l, i, j, ux, uy, vx_bc, vy_bc);
						}
					}
					write_fields2d(config, "", "buffertest", { buffer, p, p_prime }, "div, p, p_prime");
					write_x_line(0, "", "x_0", config, { buffer, p }, "div, p");

				}
				if (div_ < 1e-5 || steps > 100000)
				{
					for (unsigned int l = 0; l < N; l++)
					{
						vx[l] = ux[l];
						vy[l] = uy[l];
					}
					cout << "div = " << div_ << endl;
					//pause
					break;
				}
			}

			//solve_poisson_for_p_prime(*this, p, p0, p_bc, vx, vy, vx_bc, vy_bc);

		}
	};

	struct Projection : Configuration, Arrays
	{
		size_t iter = 0;
		walls T_bc, vx_bc, vy_bc, p_bc, C_bc, mu_bc, p_prime_bc, vx_prime_bc, vy_prime_bc, H_bc;

		Variable vx, vy, ux, uy, p, p0, T, T0, Hvisc, C, C2, C3, mu, mu2, mu3;

		bool ready = false;

		Projection() : Configuration() {}

		Projection(const Configuration& config, const Arrays& arr) : Configuration(config), Arrays(arr)
		{
			init_bc({ &T_bc, &vx_bc, &vy_bc, &p_bc, &C_bc, &mu_bc, &p_prime_bc, &vx_prime_bc, &vy_prime_bc });

			vx = Arrays::vx;
			vy = Arrays::vy;
			ux = Arrays::ux;
			uy = Arrays::uy;
			p = Arrays::p;
			p0 = Arrays::p0;
			T = Arrays::T;
			T0 = Arrays::T0;
			Hvisc = Arrays::Hvisc;

			C = Arrays::C;
			mu = Arrays::mu;

			C2 = Arrays::C2;
			mu2 = Arrays::mu2;
			C3 = Arrays::C3;
			mu3 = Arrays::mu3;

			C.set_walls(
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 });
			mu.set_walls(
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 });
			C3.wall = C2.wall = C.wall;
			mu3.wall = mu2.wall = mu.wall;

			vx.set_walls(
				{ MathBoundary::Dirichlet, 0.0 },
				{ MathBoundary::Dirichlet, 0.0 },
				{ MathBoundary::Dirichlet, 0.0 },
				{ MathBoundary::Dirichlet, 0.0 });
			vy.wall = ux.wall = uy.wall = vx.wall;


			T.set_walls(
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Dirichlet, 1.0 },
				{ MathBoundary::Dirichlet, 0.0 });
			T.set_walls(
				{ MathBoundary::Dirichlet, 0.0 },
				{ MathBoundary::Dirichlet, 1.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 });
			T0.wall = T.wall;

			//Hvisc.set_walls(
			//	{ MathBoundary::array_value, true, ny },
			//	{ MathBoundary::array_value, true, ny },
			//	{ MathBoundary::array_value, true, nx },
			//	{ MathBoundary::array_value, true, nx });

			Hvisc.set_walls(
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 },
				{ MathBoundary::Neumann, 0.0 });

			walls bottom_heating =
			{
				{ MathBoundary::Neumann, 0.0},
				{ MathBoundary::Neumann, 0.0},
				{ MathBoundary::array_value, true, nx, 1.0},
				{ MathBoundary::array_value, true, nx, 0.0}
			};

			walls side_heating =
			{
				{ MathBoundary::array_value, true, ny, 0.0 },
				{ MathBoundary::array_value, true, ny, 1.0 },
				{ MathBoundary::array_value, true, nx },
				{ MathBoundary::array_value, true, nx }
			};
			for (unsigned int i = 0; i < nx; i++)
				side_heating[south](i) = side_heating[north](i) = x_coord(*this, i);

			T0.wall = T.wall = bottom_heating;


		}

		void update_config(const Configuration& newConfig)
		{
			static_cast<Configuration&>(*this) = newConfig;
		}

		void solve_heat_box()
		{
			iter++;

			projection_convection::temperature(*this, T, T0, vx, vy);
			swap_one(*this, T, T0);

			projection_convection::quasi_velocity(*this, ux, uy, vx, vy, T);
			//quasi_velocity_fv_v2(*this, ux, uy, vx, vy, vx_bc, vy_bc); //!

			projection_convection::solve_poisson(*this, p, p0, ux, uy);
			//projection_convection::velocity_correction(*this, p, vx, vy, ux, uy, T);
			projection_convection::velocity_correction3(*this, p, vx, vy, ux, uy, T);
		}

		void solve_test()
		{
			iter++;

			projection_convection::temperature(*this, T, T0, vx, vy);
			swap_one(*this, T, T0);
			//cout << T[nx - 1 + offset * ny / 2] << endl;

			form_H_variable_viscosity(*this, Hvisc, vx, vy);

			projection_convection::quasi_velocity_variable_viscosity(*this, ux, uy, vx, vy, T, Hvisc);
			//projection_convection::quasi_velocity(*this, ux, uy, vx, vy, T);
			projection_convection::solve_poisson(*this, p, p0, ux, uy);
			projection_convection::velocity_correction3(*this, p, vx, vy, ux, uy, T);
		}

		void solve_phase_box()
		{
			iter++;

			projection_phase::quasi_velocity(*this, ux, uy, vx, vy, C, mu);
			projection_convection::solve_poisson(*this, p, p0, ux, uy);
			projection_phase::velocity_correction_phase(*this, p, vx, vy, ux, uy, C, C.wall);
			//projection_phase::velocity_correction_cut(*this, p, vx, vy, ux, uy, C);
			make_conv_buffer(*this, buffer, vx, vy, C);
		}

		void solve_phase_box_v2()
		{
			iter++;

			projection_phase::form_H_viscosity(*this, Hvisc, C);

			projection_phase::quasi_velocity_H(*this, ux, uy, vx, vy, C, mu, Hvisc);
			projection_convection::solve_poisson(*this, p, p0, ux, uy);
			projection_phase::velocity_correction_phase(*this, p, vx, vy, ux, uy, C, C.wall);
			make_conv_buffer(*this, buffer, vx, vy, C);
		}


		void solve_phase3_box()
		{
			iter++;

			projection_phase::quasi_velocity_3_phases(*this, ux, uy, vx, vy, C, mu, C2, mu2, C3, mu3);
			projection_convection::solve_poisson(*this, p, p0, ux, uy);
			projection_phase::velocity_correction_phase3(*this, p, vx, vy, ux, uy, C, C2);
			//projection_phase::velocity_correction_cut(*this, p, vx, vy, ux, uy, C);
			make_conv_buffer(*this, buffer, vx, vy, C);
			make_conv_buffer(*this, buffer2, vx, vy, C2);
		}
	};

}
















//auto dx2_border_forward = [&config](double* f, double borderF, unsigned int l)
//{
//	return (2.0 * borderF - 3.0 * f[l] + f[l + 1]) / (3.0 * pow(config.hx * 0.5, 2)); //h = half step, should be
//};
//auto dx2_border_back = [&config](double* f, double borderF, unsigned int l)
//{
//	return (2.0 * borderF - 3.0 * f[l] + f[l - 1]) / (3.0 * pow(config.hx * 0.5, 2)); //h = half step, should be
//};
//auto dy2_border_up = [&config](double* f, double borderF,  unsigned int l)
//{
//	return (2.0 * borderF - 3.0 * f[l] + f[l + config.offset]) / (3.0 * pow(config.hy * 0.5, 2)); //h = half step, should be
//};
//auto dy2_border_down = [&config](double* f, double borderF, unsigned int l)
//{
//	return (2.0 * borderF - 3.0 * f[l] + f[l - config.offset]) / (3.0 * pow(config.hy * 0.5, 2)); //h = half step, should be
//};

/*
void quasi_velocity(Configuration& config, double* ux, double* uy, double* vx, double* vy, double* f, walls& f_bc, walls& vx_bc, walls& vy_bc)
{
	double diff = 0, conv = 0;

	for (unsigned int j = 0; j < config.ny; j++) {
		for (unsigned int i = 0; i < config.nx; i++) {
			unsigned int l = i + config.offset * j;

			diff = diffusion(l, i, j, vx, vx_bc);
			conv = conv_central(l, i, j, vx, vy, vx_bc, vy_bc, vx, vx_bc);
			//conv = conv_upwind(l, i, j, vx, vy, vx_bc, vy_bc, vx, vx_bc);
			ux[l] = vx[l] + config.tau / config.dV * (diff - conv);

			diff = diffusion(l, i, j, vy, vy_bc);
			conv = conv_central(l, i, j, vx, vy, vx_bc, vy_bc, vy, vy_bc);
			//conv = conv_upwind(l, i, j, vx, vy, vx_bc, vy_bc, vy, vy_bc);
			uy[l] = vy[l] + config.tau / config.dV * (diff - conv) + config.Ra / config.Pr * config.tau * f[l];
		}
	}
}
void temperature(Configuration& config, double* f_new, double* f, double* vx, double* vy, walls& f_bc, walls& vx_bc, walls& vy_bc)
{
	double coef = 1.0 / config.Pr;
	for (unsigned int j = 0; j < config.ny; j++) {
		for (unsigned int i = 0; i < config.nx; i++) {
			unsigned int l = i + config.offset * j;
			double diff = 0.0;
			double conv = 0.0;

			//diff += (DF(east, i == config.nx - 1) - DF(west, i == 0)) * SX;
			//diff += (DF(north, j == config.ny - 1) - DF(south, j == 0)) * SY;

			//conv += (FF(east, i == config.nx - 1) * FUX(east, i == config.nx - 1) - FF(west, i == 0) * FUX(west, i == 0)) * SX;
			//conv += (FF(north, j == config.ny - 1) * FUY(north, j == config.ny - 1) - FF(south, j == 0) * FUY(south, j == 0)) * SY;

			//conv += conv_upwind(f, vx, vy, l, i, j, f_bc, vx_bc, vy_bc);

			diff = diffusion(l, i, j, f, f_bc);
			conv = conv_central(l, i, j, vx, vy, vx_bc, vy_bc, f, f_bc);
			//conv = conv_upwind(l, i, j, vx, vy, vx_bc, vy_bc, f, f_bc);

			f_new[l] = f[l] + config.tau / config.dV * (diff / config.Pr - conv);
		}
	}
}
void velocity_correction(Configuration& config, double* p, double* vx, double* vy, double* ux, double* uy, double* f, walls& T_bc)
{
	auto dx2_border_forward = [&config](double* f, double borderF, unsigned int l)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l + 1]) / (pow(config.hx * 0.5, 2)); //h = half step, should be
	};
	auto dx2_border_back = [&config](double* f, double borderF, unsigned int l)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l - 1]) / (pow(config.hx * 0.5, 2)); //h = half step, should be
	};
	auto dy2_border_up = [&config](double* f, double borderF, unsigned int l)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l + config.offset]) / (pow(config.hy * 0.5, 2)); //h = half step, should be
	};
	auto dy2_border_down = [&config](double* f, double borderF, unsigned int l)
	{
		return (borderF - 1.5 * f[l] + 0.5 * f[l - config.offset]) / (pow(config.hy * 0.5, 2)); //h = half step, should be
	};


	for (unsigned int j = 0; j < config.ny; j++) {
		for (unsigned int i = 0; i < config.nx; i++)
		{
			unsigned int l = i + config.offset * j;
			double Pe = 0, Pw = 0, Ps = 0, Pn = 0;

			if (i == 0)
			{
				if (config.domain == closed_box) Pw = p[l] - 0.5 * config.hx * dx2_border_forward(vx, 0, l);
				if (config.domain == open_tube)  Pw = config.p_in;
			}
			else
			{
				Pw = value<Side::west>(p, l);
			}

			if (i == config.nx - 1)
			{
				if (config.domain == closed_box)  Pe = p[l] + config.hx * 0.5 * dx2_border_back(vx, 0, l);
				if (config.domain == open_tube)	  Pe = config.p_out;
			}
			else
			{
				Pe = value<Side::east>(p, l);
			}


			if (j == 0)
			{
				//if (config.domain == closed_box)
				Ps = p[l] - 0.5 * config.hy * (dy2_border_up(vy, 0, l) + config.Ra / config.Pr * T_bc[south].value);
			}
			else
			{
				Ps = value<Side::south>(p, l);
			}

			if (j == config.ny - 1)
			{
				//if (config.domain == closed_box)
				Pn = p[l] + 0.5 * config.hy * (dy2_border_down(vy, 0, l) + config.Ra / config.Pr * T_bc[north].value);
			}
			else
			{
				Pn = value<Side::north>(p, l);
			}

			vx[l] = ux[l] - config.tau * (Pe - Pw) * SX / config.dV;
			vy[l] = uy[l] - config.tau * (Pn - Ps) * SY / config.dV;
		}
	}
}
*/