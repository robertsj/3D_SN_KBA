/*
 * miniapp.cc
 *
 *  Created on: Nov 4, 2014
 *      Author: kevin
 */

#include <iostream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <omp.h>
#include <string>

using namespace std;

#include "miniapp.hh"
#include "auxiliary_function.hh"

typedef std::vector<real> v_dbl;
#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

inline int Get_index(int N, int block_size, int blockID_x, int blockID_y, int blockID_z, int z_local, int x_or_y_local)
{
	return blockID_x * (N + 1) * (N + 1) * block_size * block_size + blockID_y * (N + 1) * block_size * block_size +
			blockID_z * block_size * block_size + z_local * block_size + x_or_y_local;
}

inline int Get_index(int N, int block_size, int n_eg, int blockID_x, int blockID_y, int blockID_z, int z_local, int x_or_y_local)
{
	return blockID_x * (N + 1) * (N + 1) * block_size * block_size * n_eg + blockID_y * (N + 1) * block_size * block_size * n_eg +
			blockID_z * block_size * block_size * n_eg + z_local * block_size * n_eg + x_or_y_local * n_eg;
}

Solver::Solver(int n_eg_in, int n_a_in, int cm_in, int fm_in, int upscatter_in, int iter_in)
{
	assert(upscatter_in < n_eg_in);
	n_eg = n_eg_in;
	n_a = n_a_in;
	N_A = n_a * n_a;
	cm = cm_in;
	fm = fm_in;
	upscatter = upscatter_in;
	iter = iter_in;
	cout << "input parameters \n" << "# of energy groups: " << n_eg << "\n"
			<< "# of angles in each octant: " << N_A << "\n"
			<< "# of coarse cells in each direction: " << cm << " \n"
			<< "# of fine cells in one coarse cell in each direction: "<< fm << "\n"
			<< "# of upscatter energy groups: " << upscatter << "\n"
			<< "# of running iterations: " << iter << "\n";

	totNFM_x = cm * fm;
	totNFM_y = cm * fm;
	totNFM_z = cm * fm;
	phi_size = totNFM_z * totNFM_x * totNFM_y  * n_eg;

	phi = new real [phi_size];
	Q = new real [phi_size];

	//Since the core mesh and the fine mesh is uniformly distributed over the total area, so we have the following Delta_x,Delta_y.
	//In the case the total length in X and Y direction is 1.0.
	Delta_x = 1.0 / totNFM_x;
	Delta_y = 1.0 / totNFM_y;
	Delta_z = 1.0 / totNFM_z;

	n_m = cm * cm * cm;
	//Let each part of the core mesh contain a random unique material
	srand(time(0));
	for(int i = 0; i < n_m; i++)
		RegMat.push_back(i);
	random_shuffle(RegMat.begin(),RegMat.end());

	//allocate and initiate fmmid, fine mesh (cell) material ID
	fmmid = new int [totNFM_x * totNFM_y * totNFM_z];
	for (int z = 0; z < cm; z ++)
		for(int x = 0; x < cm; x ++)
			for (int y = 0; y < cm; y ++){
				int m = RegMat[z * cm * cm + x * cm + y];
				for (int zfm = z * fm; zfm < (z + 1) * fm; zfm ++)
					for(int xfm = x * fm; xfm < (x + 1) * fm; xfm ++)
						for(int yfm = y * fm; yfm < (y + 1) * fm; yfm ++)
							fmmid[zfm * totNFM_x * totNFM_y + xfm * totNFM_y + yfm] = m;
			}


	//allocate the total cross section point and assign all to 0.5
	SigT = new real [n_m * n_eg];
	SetValue(SigT, n_m * n_eg, 0.5);
	//Allocate the scattering cross section to a low tridiagnal matrix and assign the value 0.2
	SigS = new real** [n_m];
	for(int i = 0; i < n_m; i++)
	{
		SigS[i] = new real* [n_eg];
		//for energy groups without upscattering, only energy groups with energy >=
		//(group ID <=) them can scatter into them
		for (int j = 0; j < n_eg - upscatter; j++)
			SigS[i][j] = new real [j + 1];
		//for energy groups with upscattering
		for(int j = n_eg - upscatter; j < n_eg; j++)
			SigS[i][j] = new real [n_eg];
	}

	for(int i = 0; i < n_m; i++)
	{
		for(int j = 0; j < n_eg - upscatter; j++)
			for(int k = 0; k < j + 1; k++)
				SigS[i][j][k] = 0.2;
		for(int j = n_eg - upscatter; j < n_eg; j++)
			for(int k = 0; k < n_eg; k++)
				SigS[i][j][k] = 0.2;
	}

	//Initialize mu,xi,eta;
	get_quadrature();
	cout << "--------------------------------- \n initialization finished \n--------------------------------- \n \n \n";
}

Solver::~Solver()
{
	if(phi)
		delete[] phi;

	if(Q)
		delete[] Q;

	if(SigT)
		delete[] SigT;

	if(SigS)
	{
		for(int i = 0; i < n_m; i++)
		{
			for(int j = 0; j < n_eg; j++)
				delete[] SigS[i][j];
			delete[] SigS[i];
		}
		delete[] SigS;
	}

	if(fmmid)
		delete[] fmmid;
}

void Solver::get_quadrature()
{
	/* polar angle theta; azimuthal angle phi
	 * use cosine of polar angle, omega, to represent directions, so divided [0, 1] into n_a
	 * sections;
	 * directly use azimuthal angle phi to represent directions, so divided [0, pi/2] into n_a
	 * sections */

	v_dbl phi_0(n_a, 0.0);
	real delta_phi = 0.5 * M_PI / n_a;
	for (int i = 0; i < n_a; ++i)
		phi_0[i] = real(i) * delta_phi + 0.5 * delta_phi;

	v_dbl cos_theta(n_a, 0.0);
	real delta_theta = 1.0 / n_a;
	for (int i = 0; i < n_a; ++i)
		cos_theta[i] = real(i) * delta_theta + 0.5 * delta_theta;

	for (int i = 0; i < n_a; ++i) //index of azimuthal angle
	{
		for (int j = 0; j < n_a; ++j) //index of polar angle
		{
			real sin_theta = sqrt(1 - cos_theta[j] * cos_theta[j]);
			// mu = cos(phi)*sin(theta)
			mu.push_back(cos(phi_0[i]) * sin_theta);
			// eta = sin(phi)*sin(theta)
			eta.push_back(sin(phi_0[i]) * sin_theta);
			// xi = cos(theta)
			xi.push_back(cos_theta[j]);
		}
	}
}

void Solver::Calculate(string sweepfun, int nTs_in)
{
	time_used_total = omp_get_wtime();
	nTs = nTs_in;
	cout << "# of threads is " << nTs << endl;
	N = sqrt(nTs);

	block_size = totNFM_x / N;
	int start_TID[(2 * N - 1) * 2];// start_TID -- TID info for planes having starting threads;
	Set_start_TID(start_TID, N);
	assert(sweepfun == "aes" || sweepfun == "ase" || sweepfun == "eas" || sweepfun == "esa" || sweepfun == "sae" || sweepfun == "sea");
	//Initial guess of phi, 0.0
	SetValue(phi, phi_size, 0.0);

	time_used_sweep = 0.0;

	//convergence parameters
	real eps_phi = 1.0e-5;
	real err_phi = 1.0;
	int it = 0;

	while(err_phi > eps_phi && it < iter)
	{
		real phi0[phi_size];
		copy(phi0, 0, phi, 0, phi_size);

		//Cauculate Q
		for (int k = 0; k < totNFM_z; k++)
			for (int i = 0; i < totNFM_x; i++)
			{
				for (int j = 0; j < totNFM_y; j++)
				{
					int m = fmmid[k * totNFM_x * totNFM_y + i * totNFM_y + j];
					//Here is the advantage of the allocation SigS to a low diagnal matrix
					for(int e = 0; e < n_eg - upscatter; ++e)
					{
						real totsca = 0.0;
						for(int ee = 0; ee < e + 1; ee++)
							totsca += phi[k * totNFM_x * totNFM_y * n_eg + i * totNFM_y * n_eg + j * n_eg + ee] * SigS[m][e][ee];
						//Using pi / 4 replace S[i][j][k][e] since it is uniform and
						//isotropic source.
						Q[k * totNFM_x * totNFM_y * n_eg + i * totNFM_y * n_eg + j * n_eg + e] =
								1.0 / (M_PI * 4.0) +  totsca / (M_PI * 4.0);
					}
					for(int e = n_eg - upscatter; e < n_eg; e++)
					{
						real totsca = 0.0;
						for(int ee = 0; ee < n_eg; ee++)
							totsca += phi[k * totNFM_x * totNFM_y * n_eg + i * totNFM_y * n_eg + j * n_eg + ee] * SigS[m][e][ee];
						Q[k * totNFM_x * totNFM_y * n_eg + i * totNFM_y * n_eg + j * n_eg + e] =
								1.0 / (M_PI * 4.0) +  totsca / (M_PI * 4.0);
					}
				}
			}

		real sweep_begin_time = omp_get_wtime();

		if(sweepfun == "aes")
			sweep_aes(start_TID);
		else if(sweepfun == "ase")
			sweep_ase(start_TID);
		else if(sweepfun == "eas")
			sweep_eas(start_TID);
		else if(sweepfun == "esa")
			sweep_esa(start_TID);
		else if(sweepfun == "sae")
			sweep_sae(start_TID);
		else if(sweepfun == "sea")
			sweep_sea(start_TID);

		real sweep_end_time = omp_get_wtime();
		time_used_sweep += sweep_end_time - sweep_begin_time;

		real max = fabs(phi[0] - phi0[0]) / phi0[0];
		for (int i = 1; i < phi_size; i++)
			if (fabs(phi[i] - phi0[i]) / phi0[i] > max)
				max = fabs(phi[i] - phi0[i]) / phi0[i];
		err_phi = max;
		cout << "iteration " << it << " error of phi is " << err_phi << endl;
		it = it + 1;
	}

	time_used_total = omp_get_wtime() - time_used_total;

	cout << "\n--------------------------------------\nsummary\n--------------------------------------";

	if (it <= iter && err_phi < eps_phi)
	{
		cout << "\nconverged in " << it << " iterations" << "\n";
		cout << "time used in " << sweepfun << " sweep of " << it << " iterations is " << time_used_sweep << " s" << endl;
		cout << "elapsed time per sweep is " << time_used_sweep / it << " s" << endl;
		cout << "total time is " << time_used_total << " s" << endl;

		if (phi_size < 10){
			cout << "\nthe resulting phi is \n";
			for(int i = 0; i < phi_size; i++)
				cout << phi[i] << " ";
			cout << '\n';}
		else{
			cout << "\nthe first 10 resulting phi is \n";
			for(int i = 0; i < 10; i++)
				cout << phi[i] << " ";
			cout << endl;}
		cout << "\nerr_phi is " << err_phi << "\n\n";
	}

	else
	{
		cout << "\ndo not converge and time used in "<< sweepfun << " sweep in " << iter << " iterations is " << time_used_sweep << " s" << endl;
		cout << "elapsed time per sweep is " << time_used_sweep / iter << " s" << endl;
		cout << "total time used is " << time_used_total << " s" << endl;

		if (phi_size < 10){
			cout << "\nthe resulting phi is \n";
			for (int i = 0; i < phi_size; i++)
				cout << phi[i] << " ";
			cout << endl;}
		else{
			cout << "\nthe first 10 resulting phi is \n";
			for(int i = 0; i < 10; i++)
				cout << phi[i] << " ";
			cout << endl;}
		cout << "\nerr_phi is " << err_phi << "\n\n";
	}
}


//*************************************************** KBA ***************************************************

void Solver::sweep_aes(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size];

	//Octant loop
	for(int o = 0; o < 8; o++)
	{
		DetermineDir(o, forward_z, forward_x, forward_y);
		//Angle loop
		for(int a = 0; a < N_A; a++)
		{
			const real muDelta = 2.0 * mu[a] / Delta_y;
			const real etaDelta = 2.0 * eta[a] / Delta_x;
			const real xiDelta = 2.0 * xi[a] / Delta_z;
			const real sum = muDelta + etaDelta + xiDelta;
			//Energy loop
			for(int e = 0; e < n_eg; ++e)
			{
				//for simplicity, set all to 0, in fact only the start bd_info need to set 0
				SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size, 0.0);
				SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size, 0.0);

#pragma omp parallel num_threads(nTs)
				{
					//thread-private variables
					int TID = omp_get_thread_num();
					int start_plane, end_plane;
					Find_start_plane(start_TID, N, TID, start_plane);
					end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
					//block ID
					int blockID_x, blockID_y, blockID_z;
					Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);
					real bd_info_z[block_size * block_size]; //[X][Y]
					SetValue(bd_info_z, block_size * block_size, 0.0);
					real c;

					for(int p = 0; p < 3 * N - 2; p++)
					{
						if (p >= start_plane && p < end_plane)//select working threads
						{
							blockID_z = p - start_plane;
							//block sweep
							for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
							{
								const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
								const int z_local = z_global % block_size;
								for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
								{
									const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
									const int x_local = x_global % block_size;
									const int index_x = Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, x_local);
									for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
									{
										const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
										const int y_local = y_global % block_size;
										const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
										const int index_y = Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, y_local);
										const int index_z = x_local * block_size + y_local;
										c = (muDelta * bd_info_x[index_x] +	etaDelta * bd_info_y[index_y] +	xiDelta * bd_info_z[index_z] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c;
										bd_info_x[index_x] = 2.0 * c - bd_info_x[index_x];
										bd_info_y[index_y] = 2.0 * c - bd_info_y[index_y];
										bd_info_z[index_z] = 2.0 * c - bd_info_z[index_z];
									}
								}
							}
							//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
							copy(bd_info_x, Get_index(N, block_size, blockID_x, blockID_y + 1, blockID_z, 0, 0),
									bd_info_x, Get_index(N, block_size, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size);
							copy(bd_info_y, Get_index(N, block_size, blockID_x + 1, blockID_y, blockID_z, 0, 0), bd_info_y, Get_index(N, block_size, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size);
						}
#pragma omp barrier
					}
				}
			}
		}
	}
}


void Solver::sweep_ase(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * n_eg], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * n_eg];

	//Octant loop
	for(int o = 0; o < 8; o++)
	{
		DetermineDir(o, forward_z, forward_x, forward_y);
		//Angle loop
		for(int a = 0; a < N_A; a++)
		{
			const real muDelta = 2.0 * mu[a] / Delta_y;
			const real etaDelta = 2.0 * eta[a] / Delta_x;
			const real xiDelta = 2.0 * xi[a] / Delta_z;
			const real sum = muDelta + etaDelta + xiDelta;

			//for simplicity, set all to 0, in fact only the start bd_info need to set 0
			SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * n_eg, 0.0);
			SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * n_eg, 0.0);

#pragma omp parallel num_threads(nTs)
			{
				//thread-private variables
				int TID = omp_get_thread_num();
				int start_plane, end_plane;
				Find_start_plane(start_TID, N, TID, start_plane);
				end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
				//block ID
				int blockID_x, blockID_y, blockID_z;
				Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);
				real bd_info_z[block_size * block_size * n_eg]; //[X][Y][E]
				SetValue(bd_info_z, block_size * block_size * n_eg, 0.0);
				real c[n_eg];

				for(int p = 0; p < 3 * N - 2; p++)
				{
					if (p >= start_plane && p < end_plane)//select working threads
					{
						blockID_z = p - start_plane;
						//block sweep
						for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
						{
							const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
							const int z_local = z_global % block_size;
							for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
							{
								const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
								const int x_local = x_global % block_size;
								const int index_x = Get_index(N, block_size, n_eg, blockID_x, blockID_y, blockID_z, z_local, x_local);
								for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
								{

									const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
									const int y_local = y_global % block_size;
									const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
									const int index_y = Get_index(N, block_size, n_eg, blockID_x, blockID_y, blockID_z, z_local, y_local);
									const int index_z = x_local * block_size * n_eg + y_local * n_eg;
									//Energy loop
									for(int e = 0; e < n_eg; ++e)
									{
										c[e] = (muDelta * bd_info_x[index_x + e] + etaDelta * bd_info_y[index_y + e] +	xiDelta * bd_info_z[index_z + e] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c[e];
										bd_info_x[index_x + e] = 2.0 * c[e] - bd_info_x[index_x + e];
										bd_info_y[index_y + e] = 2.0 * c[e] - bd_info_y[index_y + e];
										bd_info_z[index_z + e] = 2.0 * c[e] - bd_info_z[index_z + e];
									}
								}
							}
						}
						//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
						copy(bd_info_x, Get_index(N, block_size, n_eg, blockID_x, blockID_y + 1, blockID_z, 0, 0),
								bd_info_x, Get_index(N, block_size, n_eg, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * n_eg);
						copy(bd_info_y, Get_index(N, block_size, n_eg, blockID_x + 1, blockID_y, blockID_z, 0, 0),
								bd_info_y, Get_index(N, block_size, n_eg, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * n_eg);
					}
#pragma omp barrier
				}
			}
		}
	}
}

void Solver::sweep_eas(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size];
	//Energy loop
	for(int e = 0; e < n_eg; ++e)
	{
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			DetermineDir(o, forward_z, forward_x, forward_y);
			//Angle loop
			for(int a = 0; a < N_A; a++)
			{
				const real muDelta = 2.0 * mu[a] / Delta_y;
				const real etaDelta = 2.0 * eta[a] / Delta_x;
				const real xiDelta = 2.0 * xi[a] / Delta_z;
				const real sum = muDelta + etaDelta + xiDelta;
				//for simplicity, set all to 0, in fact only the start bd_info need to set 0
				SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size, 0.0);
				SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size, 0.0);

#pragma omp parallel num_threads(nTs)
				{
					//thread-private variables
					int TID = omp_get_thread_num();
					int start_plane, end_plane;
					Find_start_plane(start_TID, N, TID, start_plane);
					end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
					//block ID
					int blockID_x, blockID_y, blockID_z;
					Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);
					real bd_info_z[block_size * block_size]; //[X][Y]
					SetValue(bd_info_z, block_size * block_size, 0.0);
					real c;

					for(int p = 0; p < 3 * N - 2; p++)
					{
						if (p >= start_plane && p < end_plane)//select working threads
						{
							blockID_z = p - start_plane;
							//block sweep
							for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
							{
								const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
								const int z_local = z_global % block_size;
								for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
								{
									const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
									const int x_local = x_global % block_size;
									const int index_x = Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, x_local);
									for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
									{
										const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
										const int y_local = y_global % block_size;
										const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
										const int index_y = Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, y_local);
										const int index_z = x_local * block_size + y_local;
										c = (muDelta * bd_info_x[index_x] +	etaDelta * bd_info_y[index_y] +	xiDelta * bd_info_z[index_z] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c;
										bd_info_x[index_x] = 2.0 * c - bd_info_x[index_x];
										bd_info_y[index_y] = 2.0 * c - bd_info_y[index_y];
										bd_info_z[index_z] = 2.0 * c - bd_info_z[index_z];
									}
								}
							}
							//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
							copy(bd_info_x, Get_index(N, block_size, blockID_x, blockID_y + 1, blockID_z, 0, 0),
									bd_info_x, Get_index(N, block_size, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size);
							copy(bd_info_y, Get_index(N, block_size, blockID_x + 1, blockID_y, blockID_z, 0, 0), bd_info_y, Get_index(N, block_size, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size);
						}
#pragma omp barrier
					}
				}
			}
		}
	}
}

void Solver::sweep_esa(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * N_A], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * N_A];

	//Energy loop
	for(int e = 0; e < n_eg; ++e)
	{
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			DetermineDir(o, forward_z, forward_x, forward_y);
			//for simplicity, set all to 0, in fact only the start bd_info need to set 0
			SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * N_A, 0.0);
			SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * N_A, 0.0);

#pragma omp parallel num_threads(nTs)
			{
				//thread-private variables
				int TID = omp_get_thread_num();
				int start_plane, end_plane;
				Find_start_plane(start_TID, N, TID, start_plane);
				end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
				//block ID
				int blockID_x, blockID_y, blockID_z;
				Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);

				real bd_info_z[block_size * block_size * N_A]; //[X][Y][A]
				SetValue(bd_info_z, block_size * block_size * N_A, 0.0);
				real c[N_A];

				for(int p = 0; p < 3 * N - 2; p++)
				{
					if (p >= start_plane && p < end_plane)//select working threads
					{
						blockID_z = p - start_plane;
						//block sweep
						for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
						{
							const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
							const int z_local = z_global % block_size;
							for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
							{
								const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
								const int x_local = x_global % block_size;
								const int index_x = Get_index(N, block_size, N_A, blockID_x, blockID_y, blockID_z, z_local, x_local);
								for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
								{

									const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
									const int y_local = y_global % block_size;
									const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
									const int index_y = Get_index(N, block_size, N_A, blockID_x, blockID_y, blockID_z, z_local, y_local);
									const int index_z = x_local * block_size * N_A + y_local * N_A;
									//Angle loop
									for(int a = 0; a < N_A; a++)
									{
										const real muDelta = 2.0 * mu[a] / Delta_y;
										const real etaDelta = 2.0 * eta[a] / Delta_x;
										const real xiDelta = 2.0 * xi[a] / Delta_z;
										const real sum = muDelta + etaDelta + xiDelta;

										c[a] = (muDelta * bd_info_x[index_x + a] + etaDelta * bd_info_y[index_y + a] +	xiDelta * bd_info_z[index_z + a] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c[a];
										bd_info_x[index_x + a] = 2.0 * c[a] - bd_info_x[index_x + a];
										bd_info_y[index_y + a] = 2.0 * c[a] - bd_info_y[index_y + a];
										bd_info_z[index_z + a] = 2.0 * c[a] - bd_info_z[index_z + a];
									}
								}
							}
						}
						//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
						copy(bd_info_x, Get_index(N, block_size, N_A, blockID_x, blockID_y + 1, blockID_z, 0, 0),
								bd_info_x, Get_index(N, block_size, N_A, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * N_A);
						copy(bd_info_y, Get_index(N, block_size, N_A, blockID_x + 1, blockID_y, blockID_z, 0, 0),
								bd_info_y, Get_index(N, block_size, N_A, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * N_A);
					}
#pragma omp barrier
				}
			}
		}
	}
}

void Solver::sweep_sae(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	const int AE = N_A * n_eg;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * AE], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * AE];

	//Octant loop
	for(int o = 0; o < 8; o++)
	{
		DetermineDir(o, forward_z, forward_x, forward_y);
		//for simplicity, set all to 0, in fact only the start bd_info need to set 0
		SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * AE, 0.0);
		SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * AE, 0.0);

#pragma omp parallel num_threads(nTs)
		{
			//thread-private variables
			int TID = omp_get_thread_num();
			int start_plane, end_plane;
			Find_start_plane(start_TID, N, TID, start_plane);
			end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
			//block ID
			int blockID_x, blockID_y, blockID_z;
			Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);

			real bd_info_z[block_size * block_size * AE]; //[X][Y][A][E]
			SetValue(bd_info_z, block_size * block_size * AE, 0.0);
			real c[AE];


			for(int p = 0; p < 3 * N - 2; p++)
			{
				if (p >= start_plane && p < end_plane)//select working threads
				{
					blockID_z = p - start_plane;
					//block sweep
					for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
					{
						const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
						const int z_local = z_global % block_size;
						for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
						{
							const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
							const int x_local = x_global % block_size;
							const int index_x = Get_index(N, block_size, AE, blockID_x, blockID_y, blockID_z, z_local, x_local);
							for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
							{
								const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
								const int y_local = y_global % block_size;
								const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
								const int index_y = Get_index(N, block_size, AE, blockID_x, blockID_y, blockID_z, z_local, y_local);
								const int index_z = x_local * block_size * AE + y_local * AE;
								//Angle loop
								for(int a = 0; a < N_A; a++)
								{
									const real muDelta = 2.0 * mu[a] / Delta_y;
									const real etaDelta = 2.0 * eta[a] / Delta_x;
									const real xiDelta = 2.0 * xi[a] / Delta_z;
									const real sum = muDelta + etaDelta + xiDelta;
									const int A = a * n_eg;
									//Energy loop
									for(int e = 0; e < n_eg; ++e)
									{
										c[A + e] = (muDelta * bd_info_x[index_x + A + e] + etaDelta * bd_info_y[index_y + A + e] +	xiDelta * bd_info_z[index_z + A + e] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c[a * n_eg + e];
										bd_info_x[index_x + A + e] = 2.0 * c[A + e] - bd_info_x[index_x + A + e];
										bd_info_y[index_y + A + e] = 2.0 * c[A + e] - bd_info_y[index_y + A + e];
										bd_info_z[index_z + A + e] = 2.0 * c[A + e] - bd_info_z[index_z + A + e];
									}
								}
							}
						}
					}
					//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
					copy(bd_info_x, Get_index(N, block_size, AE, blockID_x, blockID_y + 1, blockID_z, 0, 0),
							bd_info_x, Get_index(N, block_size, AE, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * AE);
					copy(bd_info_y, Get_index(N, block_size, AE, blockID_x + 1, blockID_y, blockID_z, 0, 0),
							bd_info_y, Get_index(N, block_size, AE, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * AE);
				}
#pragma omp barrier
			}
		}
	}
}

void Solver::sweep_sea(int start_TID[])
{
	const real weight = 0.5 * M_PI / N_A;
	const int EA = n_eg * N_A;
	SetValue(phi, phi_size, 0.0);
	bool forward_x, forward_y, forward_z;
	//The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
	//one row is wasted along each direction
	real bd_info_x[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * EA], bd_info_y[(N + 1) * (N + 1) * (N + 1) * block_size * block_size * EA];

	//Octant loop
	for(int o = 0; o < 8; o++)
	{
		DetermineDir(o, forward_z, forward_x, forward_y);
		//for simplicity, set all to 0, in fact only the start bd_info need to set 0
		SetValue(bd_info_x, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * EA, 0.0);
		SetValue(bd_info_y, (N + 1) * (N + 1) * (N + 1) * block_size * block_size * EA, 0.0);

#pragma omp parallel num_threads(nTs)
		{
			//thread-private variables
			int TID = omp_get_thread_num();
			int start_plane, end_plane;
			Find_start_plane(start_TID, N, TID, start_plane);
			end_plane = start_plane + N; //working plane for TID is [start_plane, end_plane);
			//block ID
			int blockID_x, blockID_y, blockID_z;
			Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);

			real bd_info_z[block_size * block_size * EA]; //[X][Y][A][E]
			SetValue(bd_info_z, block_size * block_size * EA, 0.0);
			real c[EA];
			for(int p = 0; p < 3 * N - 2; p++)
			{
				if (p >= start_plane && p < end_plane)//select working threads
				{
					blockID_z = p - start_plane;
					//block sweep
					for(int z0 = blockID_z * block_size; z0 < blockID_z * block_size + block_size; z0++) //always use global ordinates
					{
						const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
						const int z_local = z_global % block_size;
						for(int x0 = blockID_x * block_size; x0 < blockID_x * block_size + block_size; x0++)
						{
							const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
							const int x_local = x_global % block_size;
							const int index_x = Get_index(N, block_size, EA, blockID_x, blockID_y, blockID_z, z_local, x_local);
							for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
							{
								const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
								const int y_local = y_global % block_size;
								const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
								const int index_y = Get_index(N, block_size, EA, blockID_x, blockID_y, blockID_z, z_local, y_local);
								const int index_z = x_local * block_size * EA + y_local * EA;
								//Energy loop
								for(int e = 0; e < n_eg; ++e)
								{
									const int E = e * N_A;
									//Angle loop
									for(int a = 0; a < N_A; a++)
									{
										const real muDelta = 2.0 * mu[a] / Delta_y;
										const real etaDelta = 2.0 * eta[a] / Delta_x;
										const real xiDelta = 2.0 * xi[a] / Delta_z;
										const real sum = muDelta + etaDelta + xiDelta;

										c[E + a] = (muDelta * bd_info_x[index_x + E + a] + etaDelta * bd_info_y[index_y + E + a] +	xiDelta * bd_info_z[index_z + E + a] +
												Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c[E + a];
										bd_info_x[index_x + E + a] = 2.0 * c[E + a] - bd_info_x[index_x + E + a];
										bd_info_y[index_y + E + a] = 2.0 * c[E + a] - bd_info_y[index_y + E + a];
										bd_info_z[index_z + E + a] = 2.0 * c[E + a] - bd_info_z[index_z + E + a];
									}
								}
							}
						}

					}
					//after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
					copy(bd_info_x, Get_index(N, block_size, EA, blockID_x, blockID_y + 1, blockID_z, 0, 0),
							bd_info_x, Get_index(N, block_size, EA, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * EA);
					copy(bd_info_y, Get_index(N, block_size, EA, blockID_x + 1, blockID_y, blockID_z, 0, 0),
							bd_info_y, Get_index(N, block_size, EA, blockID_x, blockID_y, blockID_z, 0, 0), block_size * block_size * EA);
				}
#pragma omp barrier
			}
		}
	}
}

real Solver::get_sweeptime(){
	return time_used_sweep;
}

real Solver::get_totaltime(){
	return time_used_total;
}
