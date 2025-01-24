

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
		&hostptr.omega, &hostptr.omega0, &hostptr.ksi, &hostptr.ksi0 }, host.N);
	allocate_device_arrays({ &devptr.T, &devptr.T0, &devptr.C, &devptr.C0,
		& devptr.omega, & devptr.omega0, & devptr.ksi, & devptr.ksi0 }, host.N);
	
	CuPoisson stream_poisson;
	void* args[] = { &devptr.ksi, &devptr.ksi0, &devptr.omega };
	stream_poisson.set_kernel(stream::poisson_stream, args, KERNEL2D);
	stream_poisson.set_main_field(devptr.ksi, devptr.ksi0, host.N);

	copyArrayFromHostToDevice({ devptr.T, devptr.T0 }, hostptr.T, host.N);



	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);
	check << <1, 1 >> > ();
	cudaDeviceSynchronize();
	pause




	
	deleteFilesInDirectory(L"fields");
	
	host.xbc = periodic;

reset:
	if (run.call_i >= 100) return 0;
	if (run.call_i > 0)
	{
		host.Ra += -50;
	}
	run.call_i++;

	
	double max_ksi = 0, max_omega = 0;
	Checker check_ksi(&max_ksi, &run.timeq, Checker::ExitType::Relative, "ksi");
	Checker check_omega(&max_omega, &run.timeq, Checker::ExitType::Relative, "omega");


	ftimer.start("main");

	cudaMemcpyToSymbol(dev, &host, sizeof(Configuration), 0, cudaMemcpyHostToDevice);
	while (run.stop_signal == 0)
	{
		run.iter++;
		run.timeq += host.tau;

		stream::vorticity << <KERNEL2D >> > (devptr.omega, devptr.omega0, devptr.ksi, devptr.T, devptr.C);
		stream::temperature_2d << <KERNEL2D >> > (devptr.T, devptr.T0, devptr.ksi);
		//swap_one << <KERNEL1D >> > (devptr.T0, devptr.T);
		swap_two <<< KERNEL1D>>>(devptr.omega0, devptr.omega, devptr.T0, devptr.T);

		stream_poisson.solve();

		if (run.every(10000) || run.iter == 1)
		{
			cudaMemcpy(hostptr.T, devptr.T, host.Nbytes, cudaMemcpyDeviceToHost);
			cudaMemcpy(hostptr.ksi, devptr.ksi, host.Nbytes, cudaMemcpyDeviceToHost);
			cudaMemcpy(hostptr.omega, devptr.omega, host.Nbytes, cudaMemcpyDeviceToHost);

			max_ksi = absmax(hostptr.ksi, N);
			max_omega = absmax(hostptr.omega, N);

			check_ksi.update();
			check_omega.update();


			write_fields2d(run.iter, host, { hostptr.ksi, hostptr.T }, "ksi, T");

			
			if (run.iter == 1) w_temporal << "t, time(sec), max_ksi, min_ksi" << endl;
			w_temporal << run.timeq << " " << ftimer.update_and_get("main") 
				<< " " << signedmax(hostptr.ksi, N) << " " << signedmin(hostptr.ksi, N)
				<< endl;

			cout << "t= " << run.timeq << " " << hostptr.ksi[INDEX(5, 5, 0)] << endl;
			cout << "ksi check   = " << check_ksi.dif << " " << check_ksi.dif_rel << endl;
			cout << "omega check = " << check_omega.dif << " " << check_omega.dif_rel << endl;



			if ((run.call_i == 1 && run.timeq >= 1) || (run.timeq >= 1 && run.call_i > 1))
			{
				if (check_ksi.ready_to_exit) run.stop_signal = 1;
			}
		}



		if (run.stop_signal == 1)
		{
			if (run.call_i == 1)
				w_final << "Ra, ksi_max, omega_max" << endl;
			w_final << host.Ra << " " << max_ksi << " " << max_omega
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
