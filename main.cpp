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


struct RUN_STATE
{
	int stop_signal = 0, read_database = 0, call_i = 0;
	double timeq = 0, time_sec = 0;
	int cpu_only = 0;
	size_t iter = 0;
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
	PhysicalValues phys;
	StatValues stat;
	Backup backup("recovery", true);

	std::ofstream w_final, w_temporal;

	// for the stream function - vorticity method
	allocate_host_arrays({ &hostptr.T, &hostptr.T0, &hostptr.C, &hostptr.C0, 
		&hostptr.omega, &hostptr.omega0, &hostptr.ksi, &hostptr.ksi0, &hostptr.buffer, &hostptr.buffer2, &hostptr.vx, &hostptr.vy }, host.N);
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

	if (run.read_database)
	{
		w_final.open("w_final.dat", std::ofstream::app);
		w_temporal.open("w_temporal.dat", std::ofstream::app);
		backup.read(run.iter, run.timeq, run.call_i, host, { hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C });

		copyArrayFromHostToDevice({ devptr.T, devptr.T0 }, hostptr.T, host.N);
		copyArrayFromHostToDevice({ devptr.C, devptr.C0 }, hostptr.C, host.N);
		copyArrayFromHostToDevice({ devptr.ksi, devptr.ksi0 }, hostptr.ksi, host.N);
		copyArrayFromHostToDevice({ devptr.omega, devptr.omega0 }, hostptr.omega, host.N);
	}
	else
	{
		w_final.open("w_final.dat");
		w_temporal.open("w_temporal.dat");
		deleteFilesInDirectory(L"fields");
	}

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);
	check << <1, 1 >> > ();
	cudaDeviceSynchronize();
	//pause

	if (run.read_database == 0)		stream::disturb << <1, 1 >> > (10, 10, devptr.omega0, 0.1);
	if (run.read_database == 0)		hostptr.omega0[INDEX(10,10, 0)] = 0.1;

reset:
	if (run.call_i >= 1000) return 0;
	if (run.call_i > 0 && run.iter == 0)	host.Ra += +100;
	

	Checker check_ksi(&stat.ksi_sum, &run.timeq, Checker::ExitType::Relative, "ksi", 1e-6);
	Checker check_omega(&stat.omega_sum, &run.timeq, Checker::ExitType::Relative, "omega");
	Checker check_C(&stat.C_sum, &run.timeq, Checker::ExitType::Relative, "C");
	

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);

	ftimer.start("main");
	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += host.tau;
		if (run.timeq > 100000) run.stop_signal = 2;


		ftimer.start("calc");

		if (!run.cpu_only)	{
			stream::vorticity << <KERNEL2D >> > (devptr.omega, devptr.omega0, devptr.ksi, devptr.T, devptr.C);
			stream::temperature_2d << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
			stream::concentration_2d << <KERNEL2D >> > (devptr.C, devptr.C0, devptr.ksi);
			swap_three << < KERNEL1D >> > (devptr.omega0, devptr.omega, devptr.T0, devptr.T, devptr.C0, devptr.C);
			stream_poisson.solve();
		}

		if (run.cpu_only) {
			stream_cpu::vorticity(hostptr.omega, hostptr.omega0, hostptr.ksi, hostptr.T, hostptr.C);
			stream_cpu::temperature_2d(hostptr.T, hostptr.T0, hostptr.ksi);
			stream_cpu::concentration_2d(hostptr.C, hostptr.C0, hostptr.ksi);
			stream_cpu::swap_three(hostptr.omega0, hostptr.omega, hostptr.T0, hostptr.T, hostptr.C0, hostptr.C);
			cpuPoisson.solve(hostptr.ksi, hostptr.ksi0, hostptr.omega);
		}
		ftimer.end("calc");

		// output
		//if (run.every(freq) || run.iter == 1)
		if (run.every_time(host.tau, 1) || run.iter == 1)
		{
			if (!run.cpu_only) 	{
				cudaMemcpy(hostptr.T, devptr.T, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.C, devptr.C, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.ksi, devptr.ksi, host.Nbytes, cudaMemcpyDeviceToHost);
				cudaMemcpy(hostptr.omega, devptr.omega, host.Nbytes, cudaMemcpyDeviceToHost);
			}

			stat.ksi_max = absmax(hostptr.ksi, N);
			stat.ksi_sum = sum_abs(hostptr.ksi, N);
			stat.omega_sum = sum_abs(hostptr.omega, N);
			stat.C_sum_signed = sum_signed(hostptr.C, N);
			check_ksi.update();
			check_omega.update();


			make_full(host, hostptr.buffer, hostptr.T);
			make_full(host, hostptr.buffer2, hostptr.C);
			transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);


			cout << endl << "Ra = " << host.Ra << ", t= " << run.timeq << " " << run.iter << endl;
			cout << "Stat: " << stat.ksi_sum << " " << stat.omega_sum << " " << stat.C_sum_signed << " " << hostptr.ksi[INDEX(5, 5, 0)] << endl;

			if (run.iter == 1) w_temporal << "t, time(sec), time(sec)v2, max_ksi, omega_sum, ksi_point, T_point, C_point, Csum, " << " Ra=" << host.Ra << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main") << " " << ftimer.get("calc")
				<< " " << stat.ksi_max << " " << stat.omega_sum << " " << hostptr.ksi[INDEX(10, 10, 0)] << " " << hostptr.T[INDEX(10, 10, 0)] << " " << hostptr.C[INDEX(10, 10, 0)]
				<< " " << stat.C_sum_signed 
				<< endl;



			if ((run.call_i == 0 && run.timeq >= 0) || (run.timeq >= 0 && run.call_i > 0))
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
