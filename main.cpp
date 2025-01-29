

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



struct RUN_STATE
{
	int stop_signal = 0, read_database = 0, call_i = 0;
	double timeq = 0, time_sec = 0;
	size_t iter = 0;
	bool every(int n) { return iter % n == 0; }
	bool every_time(double tau, double t) {int tt = (int)(round(1.0 / tau) * t); return iter % (tt) == 0;}
	void reset()
	{
		timeq = 0;
		iter = 0;
	}
} run;

__global__ void check2(double* f, unsigned int l)
{
	printf("check: %f \n", f[l]);
}

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
	Backup backup;

	std::ofstream w_final, w_temporal;


	if (run.read_database)
	{
		w_final.open("w_final.dat", std::ofstream::app);
		w_temporal.open("w_temporal.dat", std::ofstream::app);
	}
	else
	{
		w_final.open("w_final.dat");
		w_temporal.open("w_temporal.dat");
	}

	// for the stream function - vorticity method
	allocate_host_arrays({ &hostptr.T, &hostptr.T0, &hostptr.C, &hostptr.C0, 
		&hostptr.omega, &hostptr.omega0, &hostptr.ksi, &hostptr.ksi0, &hostptr.buffer, &hostptr.buffer2, &hostptr.vx, &hostptr.vy }, host.N);
	allocate_device_arrays({ &devptr.T, &devptr.T0, &devptr.C, &devptr.C0,
		& devptr.omega, & devptr.omega0, & devptr.ksi, & devptr.ksi0 }, host.N);
	
	CuPoisson stream_poisson;
	void* args[] = { &devptr.ksi, &devptr.ksi0, &devptr.omega };
	stream_poisson.set_kernel(stream::poisson_stream, args, KERNEL2D);
	stream_poisson.set_main_field(devptr.ksi, devptr.ksi0, host.N);

	copyArrayFromHostToDevice({ devptr.T, devptr.T0 }, hostptr.T, host.N);








	
	deleteFilesInDirectory(L"fields");
	host.xbc = periodic;

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);
	check << <1, 1 >> > ();
	cudaDeviceSynchronize();
	//pause





	stream::disturb << <1, 1 >> > (5, 5, devptr.omega0, 1);
reset:
	if (run.call_i >= 100) return 0;
	if (run.call_i > 0)
	{
		host.Ra += -5;
	}
	run.call_i++;

	
	
	Checker check_ksi(&stat.ksi_sum, &run.timeq, Checker::ExitType::Relative, "ksi", 1e-7);
	Checker check_omega(&stat.omega_sum, &run.timeq, Checker::ExitType::Relative, "omega");
	Checker check_C(&stat.C_sum, &run.timeq, Checker::ExitType::Relative, "C");
	

	ftimer.start("main");

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);

	//stream::disturb << <1, 1 >> > (5, 5, devptr.omega0, 1);
	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += host.tau;

		stream::vorticity << <KERNEL2D >> > (devptr.omega, devptr.omega0, devptr.ksi, devptr.T, devptr.C);
		stream::temperature_2d << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
		stream::concentration_2d << <KERNEL2D >> > (devptr.C, devptr.C0, devptr.ksi);

		//swap_one << <KERNEL1D >> > (devptr.T0, devptr.T);
		swap_three << < KERNEL1D >> > (devptr.omega0, devptr.omega, devptr.T0, devptr.T, devptr.C0, devptr.C);

		stream_poisson.solve();



		// output
		unsigned int freq = 1;
		//if (run.every(freq) || run.iter == 1)
		if (run.every_time(host.tau, freq) || run.iter == 1)
		{
			cudaMemcpy(hostptr.T, devptr.T, host.Nbytes, cudaMemcpyDeviceToHost);
			cudaMemcpy(hostptr.C, devptr.C, host.Nbytes, cudaMemcpyDeviceToHost);
			cudaMemcpy(hostptr.ksi, devptr.ksi, host.Nbytes, cudaMemcpyDeviceToHost);
			cudaMemcpy(hostptr.omega, devptr.omega, host.Nbytes, cudaMemcpyDeviceToHost);

			stat.ksi_max = absmax(hostptr.ksi, N);
			stat.ksi_sum = sum_abs(hostptr.ksi, N);
			stat.omega_sum = sum_abs(hostptr.omega, N);
			stat.C_sum_signed = sum_signed(hostptr.C, N);
			check_ksi.update();
			check_omega.update();


			make_full(host, hostptr.buffer, hostptr.T);
			make_full(host, hostptr.buffer2, hostptr.C);
			transform_to_velocity(host, hostptr.ksi, hostptr.vx, hostptr.vy);


			cout << endl << "Ra = " << host.Ra << ", t= " << run.timeq << endl;
			cout << "Stat: " << stat.ksi_sum << " " << stat.omega_sum << " " << stat.C_sum_signed << " " << hostptr.ksi[INDEX(5, 5, 0)] << endl;



			if (run.iter == 1) w_temporal << "t, time(sec), max_ksi, omega_sum, ksi_point, Csum, " << " Ra=" << host.Ra << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main")
				<< " " << stat.ksi_max << " " << stat.omega_sum << " " << hostptr.ksi[INDEX(10, 10, 0)]
				<< " " << stat.C_sum_signed 
				<< endl;


			if ((run.call_i == 1 && run.timeq >= 1) || (run.timeq >= 1 && run.call_i > 1))
			{
				cout << "ksi check   = " << check_ksi.dif << " " << check_ksi.dif_rel << endl;
				cout << "omega check = " << check_omega.dif << " " << check_omega.dif_rel << endl;
				if (check_ksi.ready_to_exit) run.stop_signal = 1;
			}

			if (run.every_time(host.tau, freq * 10))
			{
				write_fields2d(run.iter, host,
					{ hostptr.ksi, hostptr.T, hostptr.C, hostptr.buffer, hostptr.buffer2, hostptr.vx, hostptr.vy },
					"ksi, T, C, Tfull, Cfull, vx, vy");
			}
			if (run.every_time(host.tau, freq * 50))
			{
				backup.save(run.iter, run.timeq, host,
					{ hostptr.ksi, hostptr.omega, hostptr.T, hostptr.C }, 
					"ksi, omega, T, C");
			}
		}



		if (run.stop_signal == 1)
		{
			if (run.call_i == 1)	w_final << "Ra, ksi_max, ksi_max2, omega_sum, time(sec), t" << endl;

			w_final << host.Ra << " " << stat.ksi_max << " " << pow(stat.ksi_max, 2) << " " << stat.omega_sum
				<< " " << ftimer.update_and_get("main") << " " << run.timeq
				<< endl;

			//back.write_fields({ "vx", "vy", "vz", "C", "T", "p" }, { vx_h, vy_h, vz_h, C2_h, T_h, p_h }, Nbytes);
			//back.write_parameters(iter, timeq, host.Ra);
			//pause
			run.iter = 0;
			run.timeq = 0;
			run.stop_signal = 0;
			goto reset;
		}
		

	}
	return 0;
}
