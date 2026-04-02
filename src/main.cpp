

#include "Extras.h"
#include "ExtrasCuda.h"

#include "types_project.h"
#include "kernels/stream.h"
//#include "kernels/fv_projection.h"

#include "CuPoisson/CuPoisson.h"
#include "cuda_runtime.h"
//#include <device_launch_parameters.h>
#include <vector>
#include <iostream>


using std::cout;
using std::endl;

__constant__ Configuration dev;
Configuration host;
#include "init.h"
#include "CPUSolvers/cpu_stream.h"
#include "CPUSolvers/cpu_stream_impl.h"
#include "kernels/projection/Header.h"
//#include "CG.h"

//#include "CuBiCGStab/CuCG.h"


RUN_STATE run;

int main(int argc, char** argv)
{
	ReadingFile par("parameters.txt");
	par.reading<int>(run.read_recovery, "continue", 0); // 0 - full start, 1 - full continue, 2 - read fields, state is clear
	par.reading<int>(run.compute_mode, "cpu", 1);
	par.reading<double>(run.timeq_limit, "time_limit", 20000);
	par.reading<double>(run.timeq_minimal, "time_minimal", 0);
	par.reading<int>(run.call_max, "call_max", 100);
	par.reading<int>(run.full_reset, "full_reset", 0);
	par.reading_string(run.note, "note", "");

	GPU_ gpu(0);
	init_parameters(host);
	unsigned int &N = host.N;

	Arrays hostptr, devptr;
	CudaLaunchSetup launch(host.N, host.nx, host.ny, host.nz);
	#define KERNEL1D launch.grid1d, launch.block1d
	#define KERNEL2D launch.grid2d, launch.block2d
	FuncTimer ftimer;
	//PhysicalValues phys;
	StatValues stat;
	Backup backup("recovery", true);
	std::ofstream w_final, w_temporal;

	// for the stream function - vorticity method
	allocate_host_arrays({ &hostptr.T, &hostptr.T0, &hostptr.C, &hostptr.C0, 
		&hostptr.omega, &hostptr.omega0, &hostptr.ksi, &hostptr.ksi0, &hostptr.buffer, &hostptr.buffer2, &hostptr.vx, &hostptr.vy, &hostptr.rhs }, host.N);
	allocate_device_arrays({ &devptr.T, &devptr.T0, &devptr.C, &devptr.C0,
		& devptr.omega, & devptr.omega0, & devptr.ksi, & devptr.ksi0 }, host.N);
	
	CuPoisson stream_poisson;
	void* args[] = { &devptr.ksi, &devptr.ksi0, &devptr.omega };
	CuPoisson poisson(N, stream_cuda::poisson_stream, args, KERNEL2D);
	//stream_poisson.set_kernel(stream_cuda::poisson_stream, args, KERNEL2D);
	//stream_poisson.set_main_field(devptr.ksi, devptr.ksi0, host.N);

	stream_cpu::CuPoisson cpuPoisson;
	stream_cpu::ImplicitStream IS(N, host);

	if (run.read_recovery)
	{
		int error = backup.read(run.iter, run.timeq, run.call_i, host, { hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C }, run.read_recovery == 1);

		if (error == 0)
		{
			copyArrayFromHostToDevice({ devptr.T, devptr.T0 }, hostptr.T, host.N);
			copyArrayFromHostToDevice({ devptr.C, devptr.C0 }, hostptr.C, host.N);
			copyArrayFromHostToDevice({ devptr.ksi, devptr.ksi0 }, hostptr.ksi, host.N);
			copyArrayFromHostToDevice({ devptr.omega, devptr.omega0 }, hostptr.omega, host.N);
			if (run.compute_mode > 0)
			{
				for (unsigned int l = 0; l < host.N; l++)
				{
					hostptr.ksi0[l] = hostptr.ksi[l];
					hostptr.omega0[l] = hostptr.omega[l];
					hostptr.T0[l] = hostptr.T[l];
					hostptr.C0[l] = hostptr.C[l];
				}
			}
			auto open_type = run.read_recovery == 1 ? std::ofstream::app : std::ofstream::out;
			w_final.open("w_final.dat", open_type);
			w_temporal.open("w_temporal.dat", open_type);
			if (run.read_recovery == 2) deleteFilesInDirectory(L"fields");
		}
		if (error == 1)
		{
			run.read_recovery = 0;
		}
	}
	if (!run.read_recovery)
	{
		w_final.open("w_final.dat");
		w_temporal.open("w_temporal.dat");
		deleteFilesInDirectory(L"fields");
		init_fields(host, hostptr, devptr);
	}

	Trajectory tr(host, 1, run.read_recovery == 1);
	if (run.read_recovery == 0) { tr.x[0] = 0.2; tr.y[0] = 0.2; tr.z[0] = 0; }

reset:
	if (run.call_i >= run.call_max) return 0;
	if (run.call_i > 0 && run.iter == 0)	host.Ra += host.incr_parameter;
	if (run.full_reset == 1 && run.iter == 0)  init_fields(host, hostptr, devptr);


	Checker check_ksi(&stat.ksi_sum, &run.timeq, Checker::ExitType::Relative, "ksi", 1e-6);
	Checker check_omega(&stat.omega_sum, &run.timeq, Checker::ExitType::Relative, "omega");
	Checker check_C(&stat.C_sum, &run.timeq, Checker::ExitType::Relative, "C");

	Checker check_ksi_max(&stat.ksi_max, &run.timeq, Checker::ExitType::NoExit, "ksimax");
	Checker check_Vmax(&stat.Vmax, &run.timeq, Checker::ExitType::NoExit, "check_Vmax");
	Checker check_Ek(&stat.Ek, &run.timeq, Checker::ExitType::NoExit, "check_Ek");
	

	if (host.Ra == 0) host.Ra += host.incr_parameter;
	host.K = host.Rad / host.Ra / host.Le;
	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);
	//check << <1, 1 >> > ();	cudaDeviceSynchronize(); //pause
	
	ftimer.start("main");

	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += host.tau;
		if (run.timeq > run.timeq_limit) run.stop_signal = 2;


		ftimer.start("calc");

		if (run.compute_mode == 0)	{
			stream_cuda::vorticity << <KERNEL2D >> > (devptr.omega, devptr.omega0, devptr.ksi, devptr.T, devptr.C);
			//stream_cuda::temperature_2d << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
			//stream_cuda::temperature_2d_flux << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
			//stream_cuda::concentration_2d << <KERNEL2D >> > (devptr.C, devptr.C0, devptr.ksi);
			stream_cuda::temperature_2d_flux_full << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
			//stream_cuda::concentration_2d_full << <KERNEL2D >> > (devptr.C, devptr.C0, devptr.ksi);
			swap_three << < KERNEL1D >> > (devptr.omega0, devptr.omega, devptr.T0, devptr.T, devptr.C0, devptr.C);
			poisson.solve();
		}

		if (run.compute_mode == 1) // explicit
		{
			stream_cpu::vorticity(hostptr.omega, hostptr.omega0, hostptr.ksi, hostptr.T, hostptr.C);
			stream_cpu::temperature_2d(hostptr.T, hostptr.T0, hostptr.ksi);
			stream_cpu::concentration_2d(hostptr.C, hostptr.C0, hostptr.ksi);
			//stream_cpu::temperature_2d_full(hostptr.T, hostptr.T0, hostptr.ksi);
			//stream_cpu::temperature_2d_flux(hostptr.T, hostptr.T0, hostptr.ksi);
			//stream_cpu::concentration_2d_full(hostptr.C, hostptr.C0, hostptr.ksi);
			
			//stream_cpu::vorticity_Soret(hostptr.omega, hostptr.omega0, hostptr.ksi, hostptr.T, hostptr.C);
			//stream_cpu::temperature_2d_flux_full(hostptr.T, hostptr.T0, hostptr.ksi);
			//stream_cpu::concentration_2d_full_Soret(hostptr.C, hostptr.C0, hostptr.T0, hostptr.ksi);

			stream_cpu::swap_three(hostptr.omega0, hostptr.omega, hostptr.T0, hostptr.T, hostptr.C0, hostptr.C);

			cpuPoisson.solve(hostptr.ksi, hostptr.ksi0, hostptr.omega);
		}
		if (run.compute_mode == 2) //implicit
		{
			stream_cpu::transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);

			stream_cpu::form_rhs_vorticity(host, hostptr.rhs, hostptr.T, hostptr.C, hostptr.ksi, hostptr.omega0, hostptr.vx, hostptr.vy, true, host.tau);
			stream_cpu::solveJacobi(hostptr.omega, hostptr.omega0, hostptr.rhs, N, IS.SMO);

			stream_cpu::form_rhs_temperature(host, hostptr.rhs, hostptr.T, hostptr.vx, hostptr.vy, true, host.tau);
			stream_cpu::solveJacobi(hostptr.T, hostptr.T0, hostptr.rhs, N, IS.SMT);

			stream_cpu::form_rhs_concentration(host, hostptr.rhs, hostptr.C, hostptr.vx, hostptr.vy, true, host.tau);
			stream_cpu::solveJacobi(hostptr.C, hostptr.C0, hostptr.rhs, N, IS.SMC);

			stream_cpu::form_rhs_stream(host, hostptr.rhs, hostptr.omega, true);
			stream_cpu::solveJacobi(hostptr.ksi, hostptr.ksi0, hostptr.rhs, N, IS.SMP);
		}


		ftimer.end("calc");
		// OUTPUT
		tr.trace_all(hostptr.vx, hostptr.vy, nullptr);
		if (run.iter == 1 
			|| (run.every_time(host.tau, 1)  && run.timeq < 1000) 
			|| (run.every_time(host.tau, 10.0) && run.timeq >= 1000))
		{
			if (run.compute_mode == 0) 	{
				cudaMemcpy(hostptr.T, devptr.T, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.C, devptr.C, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.ksi, devptr.ksi, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.omega, devptr.omega, host.Nbytes, cudaMemcpyDeviceToHost);
			}

			//Nu_y(host, hostptr.T, stat.NuTop, stat.NuDown, host.heatflux ? 1.0 : 0);
			//Nu_y(host, hostptr.C, stat.ShrTop, stat.ShrDown, host.heatflux ? 1.0 : 0);





			make_full(host, hostptr.buffer, hostptr.T);
			make_full(host, hostptr.buffer2, hostptr.C);
			stream_cpu::transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);
			velocity_stats(host, hostptr.vx, hostptr.vy, stat);
			Nu_y_for_fixed_flux(host, hostptr.buffer, stat.Nu);
			Nu_y_for_fixed_flux(host, hostptr.buffer2, stat.Shr);

			stat.ksi_max = absmax(hostptr.ksi, N);
			stat.ksi_sum = sum_abs(hostptr.ksi, N);
			stat.omega_sum = sum_abs(hostptr.omega, N);
			stat.C_sum_signed = sum_signed(hostptr.C, N);
			check_ksi.update();
			check_omega.update();
			check_ksi_max.update();
			check_Vmax.update();
			check_Ek.update();

			stat.ksi_incr = check_ksi_max.ln_time; //check it


			//if (run.every_time(host.tau, 0.1))
			{
				cout << endl << "Ra = " << host.Ra << ", t= " << run.timeq << ", " << run.iter << endl;
				cout << "ksi= " << stat.ksi_max << ", Vmax= " << stat.Vmax << ", dC= " << hostptr.C[INDEX(host.nx / 2, host.ny, 0)] - hostptr.C[INDEX(host.nx / 2, 0, 0)] << ", Pe = " << stat.Pe << ", Csum = " << stat.C_sum_signed << ", T_fix = " << hostptr.T[INDEX(5, 5, 0)] << endl;
				if (!run.note.empty()) cout << "note: " << run.note << endl;
			}

			if (run.iter == 1) w_temporal << "t, time(sec), time(sec)v2, max_ksi, omega_sum, ksi_point, T_point, C_point, Vmax, Ek, Csum, dC, dC1, Nu, Shr, ksi_incr, check_Vmax, Ek_incr" << " Ra=" << host.Ra << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main") << " " << ftimer.get("calc")
				<< " " << stat.ksi_max << " " << stat.omega_sum << " " << hostptr.ksi[INDEX(10, 10, 0)] << " " << hostptr.T[INDEX(10, 10, 0)] << " " << hostptr.C[INDEX(10, 10, 0)]
				<< " " << stat.Vmax << " " << stat.Ek
				<< " " << stat.C_sum_signed 
				<< " " << hostptr.C[INDEX(host.nx / 2, host.ny, 0)] - hostptr.C[INDEX(host.nx / 2, 0, 0)]
				<< " " << hostptr.C[INDEX(host.nx / 2, host.ny - 1, 0)] - hostptr.C[INDEX(host.nx / 2, 1, 0)]
				<< " " << stat.Nu << " " << stat.Shr << " "
				<< " " << check_ksi_max.ln_time << " " << check_Vmax.ln_time << " " << check_Ek.ln_time
				<< endl;

			tr.write(run.timeq);

			if (run.timeq >= run.timeq_minimal)
			{
				cout << "ksi check   = " << check_ksi.dif << " " << check_ksi.dif_rel << endl;
				cout << "omega check = " << check_omega.dif << " " << check_omega.dif_rel << endl;
				if (check_ksi.ready_to_exit) run.stop_signal = 1;
			}

			if ((run.every_time(host.tau, 20.0) && run.timeq < 1000)
				|| (run.every_time(host.tau, 100.0) && run.timeq >= 1000))
			{
				write_fields2d(_path("fields"), _str(host.Ra) + " " + _str(run.timeq), host,
					{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C, hostptr.buffer, hostptr.buffer2, hostptr.vx, hostptr.vy },
					"ksi, omega, T, C, Tfull, Cfull, vx, vy"
				);
			}
			if (run.every_time(host.tau, 20))
			{
				backup.save(run.iter, run.timeq, run.call_i, host, { hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C }, "ksi, omega, T, C");
			}
		}

		if (run.stop_signal > 0)
		{
			if (run.call_i == 0)	w_final << "Ra, ksi_max, ksi_max2, omega_sum, Vmax, Ek, C_sum, time(sec), t, Nu, Shr" << endl;

			w_final << host.Ra << " " << stat.ksi_max << " " << pow(stat.ksi_max, 2) << " " << stat.omega_sum << " " << stat.Vmax << " " << stat.Ek << " " << stat.C_sum_signed
				<< " " << ftimer.update_and_get("main") << " " << run.timeq
				<< " " << stat.Nu << " " << stat.Shr
				<< endl;
			backup.save(run.iter, run.timeq, run.call_i, host,	{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C },	"ksi, omega, T, C");
			
			write_fields2d(_path("final"), _str(host.Ra) + " " + _str(run.timeq), host,
				{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C, hostptr.buffer, hostptr.buffer2, hostptr.vx, hostptr.vy },
				"ksi, omega, T, C, Tfull, Cfull, vx, vy");

			if (run.stop_signal == 1)
			{
				run.iter = 0; run.timeq = 0; run.stop_signal = 0; 
				run.call_i++;
				goto reset;
			}
			if (run.stop_signal == 2) break;
		}

		if (run.stop_signal == -1)
		{
			break;
		}
	} cout << "run.stop: " << run.stop_signal << endl;

	

	
	//std::ofstream w("test.dat", std::ofstream::app);
	//cout << ftimer.get("calc") << " " << host.N << " " << host.Lx << " " << run.timeq << endl;
	//w << ftimer.get("calc") << " " << host.N << " " << host.Lx << " " << run.timeq << endl;

	return 0;
}
