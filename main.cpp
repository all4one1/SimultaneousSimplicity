#include "Extras.h"
#include "ExtrasCuda.h"
#include "types_project.h"
#include "kernels/stream.h"
#include "CuPoisson/CuPoisson.h"
#include "cuda_runtime.h"
//#include <device_launch_parameters.h>
#include <vector>
#include <iostream>

using std::cout;
using std::endl;

__constant__ Configuration dev;
Configuration host;
#include "cpu_stream.h"
#include "cpu_stream_impl.h"
#include "cpusolver.h"
#include "CG.h"

struct RUN_STATE
{
	int stop_signal = 0, read_database = 0, call_i = 0, cpu_only = 0;
	double timeq = 0, time_sec = 0, timeq_limit = 10000, timeq_minimal = 0;
	size_t iter = 0;
	std::string note = "";
	bool every(int n) { return iter % n == 0; }
	bool every_time(double tau, double t) {int tt = (int)(round(1.0 / tau) * t); return iter % (tt) == 0;}
	void reset()
	{
		iter = 0;
		timeq = 0;
		stop_signal = 0;
	}
} run;

int main(int argc, char** argv)
{
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
	stream_poisson.set_kernel(stream::poisson_stream, args, KERNEL2D);
	stream_poisson.set_main_field(devptr.ksi, devptr.ksi0, host.N);

	stream_cpu::CuPoisson cpuPoisson;
	
	ReadingFile par("parameters.txt");	
	par.reading<int>(run.read_database, "continue", 0);
	par.reading<int>(run.cpu_only, "cpu", 0);
	par.reading<double>(run.timeq_limit, "time_limit", 20000);
	par.reading<double>(run.timeq_minimal, "time_minimal", 0);
	par.reading_string(run.note, "note", "");

	Trajectory tr(host, 1, run.read_database);
	tr.x[0] = 0.25;
	tr.y[0] = 0.25;
	tr.z[0] = 0;


	if (run.read_database)
	{
		w_final.open("w_final.dat", std::ofstream::app);
		w_temporal.open("w_temporal.dat", std::ofstream::app);
		backup.read(run.iter, run.timeq, run.call_i, host, { hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C });
		
		copyArrayFromHostToDevice({ devptr.T, devptr.T0 }, hostptr.T, host.N);
		copyArrayFromHostToDevice({ devptr.C, devptr.C0 }, hostptr.C, host.N);
		copyArrayFromHostToDevice({ devptr.ksi, devptr.ksi0 }, hostptr.ksi, host.N);
		copyArrayFromHostToDevice({ devptr.omega, devptr.omega0 }, hostptr.omega, host.N);

		if (run.cpu_only)
		{
			for (unsigned int l = 0; l < host.N; l++)
			{
				hostptr.ksi0[l] = hostptr.ksi[l];
				hostptr.omega0[l] = hostptr.omega[l];
				hostptr.T0[l] = hostptr.T[l];
				hostptr.C0[l] = hostptr.C[l];
 			}
		}
	}
	else
	{
		w_final.open("w_final.dat");
		w_temporal.open("w_temporal.dat");
		deleteFilesInDirectory(L"fields");
	}

	//cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);	check << <1, 1 >> > ();	cudaDeviceSynchronize();

	if (run.read_database == 0)		stream::disturb << <1, 1 >> > (10, 10, devptr.omega0, 0.1);
	if (run.read_database == 0)		hostptr.omega0[INDEX(10,10, 0)] = 0.1;

reset:
	if (run.call_i >= 1) return 0;
	if (run.call_i > 0 && run.iter == 0)	host.Ra += -100;
	

	Checker check_ksi(&stat.ksi_sum, &run.timeq, Checker::ExitType::Relative, "ksi", 1e-5);
	Checker check_omega(&stat.omega_sum, &run.timeq, Checker::ExitType::Relative, "omega");
	Checker check_C(&stat.C_sum, &run.timeq, Checker::ExitType::Relative, "C");
	

	IterativeSolver isolv;
	SparseMatrix SM(N);
	stream_cpu::ImplicitStream IS(N, host);
	BICGSTAB bi(N);



	//form_matrix_finite_difference(SM, host);
	form_matrix_finite_difference_v2(SM, host);

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);

	ftimer.start("main");
	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += host.tau;
		if (run.timeq > run.timeq_limit) run.stop_signal = 2;


		ftimer.start("calc");

		if (!run.cpu_only)	{
			stream::vorticity << <KERNEL2D >> > (devptr.omega, devptr.omega0, devptr.ksi, devptr.T, devptr.C);
			stream::temperature_2d << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
			stream::concentration_2d << <KERNEL2D >> > (devptr.C, devptr.C0, devptr.ksi);
			swap_three << < KERNEL1D >> > (devptr.omega0, devptr.omega, devptr.T0, devptr.T, devptr.C0, devptr.C);
			stream_poisson.solve();
		}

		if (run.cpu_only == 1) {
			stream_cpu::vorticity(hostptr.omega, hostptr.omega0, hostptr.ksi, hostptr.T, hostptr.C);
			stream_cpu::temperature_2d(hostptr.T, hostptr.T0, hostptr.ksi);
			stream_cpu::concentration_2d(hostptr.C, hostptr.C0, hostptr.ksi);
			stream_cpu::swap_three(hostptr.omega0, hostptr.omega, hostptr.T0, hostptr.T, hostptr.C0, hostptr.C);
			cpuPoisson.solve(hostptr.ksi, hostptr.ksi0, hostptr.omega);
		}
		if (run.cpu_only == 2)
		{
			stream_cpu::transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);

			stream_cpu::form_rhs_vorticity(host, hostptr.rhs, hostptr.T, hostptr.C, hostptr.ksi, hostptr.omega0, hostptr.vx, hostptr.vy, true, host.tau);
			stream_cpu::solveJacobi(hostptr.omega, hostptr.omega0, hostptr.rhs, N, IS.SMO);

			stream_cpu::form_rhs_temperature(host, hostptr.rhs, hostptr.T0, hostptr.vx, hostptr.vy, true, host.tau);
			stream_cpu::solveJacobi(hostptr.T, hostptr.T0, hostptr.rhs, N, IS.SMT);

			//stream_cpu::temperature_2d(hostptr.T, hostptr.T0, hostptr.ksi);
			//stream_cpu::concentration_2d(hostptr.C, hostptr.C0, hostptr.ksi);
			stream_cpu::swap_three(hostptr.omega0, hostptr.omega, hostptr.T0, hostptr.T, hostptr.C0, hostptr.C);

			stream_cpu::form_rhs_stream(host, hostptr.rhs, hostptr.omega, true);
			stream_cpu::solveJacobi(hostptr.ksi, hostptr.ksi0, hostptr.rhs, N, IS.SMP);
		}



		ftimer.end("calc");

		// OUTPUT
		if (run.every_time(host.tau, 1) || run.iter == 1)
		{
			if (!run.cpu_only) 	{
				cudaMemcpy(hostptr.T, devptr.T, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.C, devptr.C, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.ksi, devptr.ksi, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.omega, devptr.omega, host.Nbytes, cudaMemcpyDeviceToHost);
			}

			Nu_y(host, hostptr.T, stat.NuTop, stat.NuDown);
			Nu_y(host, hostptr.T, stat.ShrTop, stat.ShrDown);

			stat.ksi_max = absmax(hostptr.ksi, N);
			stat.ksi_sum = sum_abs(hostptr.ksi, N);
			stat.omega_sum = sum_abs(hostptr.omega, N);
			stat.C_sum_signed = sum_signed(hostptr.C, N);
			check_ksi.update();
			check_omega.update();

			make_full(host, hostptr.buffer, hostptr.T);
			make_full(host, hostptr.buffer2, hostptr.C);
			stream_cpu::transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);
			velocity_stats(host, hostptr.vx, hostptr.vy, stat);

			cout << endl << "Ra = " << host.Ra << ", t= " << run.timeq << ", " << run.iter << endl;
			cout << "ksi= " << stat.ksi_max << ", Vmax= " << stat.Vmax << ", Cu = " << stat.Cu << ", Pe = " << stat.Pe << ", Csum = " << stat.C_sum_signed << ", ksi_fix = " << hostptr.ksi[INDEX(5, 5, 0)] << endl;
			if (!run.note.empty()) cout << "note: " << run.note << endl;

			if (run.iter == 1) w_temporal << "t, time(sec), time(sec)v2, max_ksi, omega_sum, ksi_point, T_point, C_point, Csum, NuTop, NuDown, ShrTop, ShrDown, " << " Ra=" << host.Ra << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main") << " " << ftimer.get("calc")
				<< " " << stat.ksi_max << " " << stat.omega_sum << " " << hostptr.ksi[INDEX(10, 10, 0)] << " " << hostptr.T[INDEX(10, 10, 0)] << " " << hostptr.C[INDEX(10, 10, 0)]
				<< " " << stat.C_sum_signed 
				<< " " << stat.NuTop << " " << stat.NuDown << " " << stat.ShrTop << " " << stat.ShrDown
				<< endl;

			tr.trace_all(run.timeq, hostptr.vx, hostptr.vy, nullptr);


			if (run.timeq >= run.timeq_minimal)
			{
				cout << "ksi check   = " << check_ksi.dif << " " << check_ksi.dif_rel << endl;
				cout << "omega check = " << check_omega.dif << " " << check_omega.dif_rel << endl;
				if (check_ksi.ready_to_exit) run.stop_signal = 1;
			}

			if (run.every_time(host.tau, 20))
			{
				write_fields2d(run.iter, host,
					{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C, hostptr.buffer, hostptr.buffer2, hostptr.vx, hostptr.vy },
					"ksi, omega, T, C, Tfull, Cfull, vx, vy");
			}
			if (run.every_time(host.tau, 20))
			{
				backup.save(run.iter, run.timeq, run.call_i, host, { hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C }, "ksi, omega, T, C");
			}
		}

		if (run.stop_signal > 0)
		{
			if (run.call_i == 0)	w_final << "Ra, ksi_max, ksi_max2, omega_sum, C_sum, time(sec), t" << endl;

			w_final << host.Ra << " " << stat.ksi_max << " " << pow(stat.ksi_max, 2) << " " << stat.omega_sum << " " << stat.C_sum_signed
				<< " " << ftimer.update_and_get("main") << " " << run.timeq
				<< endl;
			backup.save(run.iter, run.timeq, run.call_i, host,	{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C },	"ksi, omega, T, C");
			
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
		

	}
	return 0;
}
