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

Solver::Solver(int n_eg_in, int n_a_in, int cm_in, int fm_in, int upscatter_in, int iter_in)
{
	assert(upscatter_in < n_eg_in);
	n_eg = n_eg_in;
	n_a = n_a_in;
	cm = cm_in;
	fm = fm_in;
	upscatter = upscatter_in;
	iter = iter_in;
	cout << "input parameters \n" << "# of energy groups: " << n_eg << "\n"
			<< "# of angles in each octant: " << n_a * n_a << "\n"
			<< "# of coarse cells in each direction: " << cm << " \n"
			<< "# of fine cells in one coarse cell in each direction: "<< fm << "\n"
			<< "# of upscatter energy groups: " << upscatter << "\n"
			<< "# of running iterations: " << iter << "\n";

	N_A = n_a * n_a;
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
	nTs = nTs_in;
	cout << "# of threads is " << nTs << endl;
	N = sqrt(nTs);

	block_size = totNFM_x / N;
	int start_TID[(2 * N - 1) * 2];// start_TID -- TID info for planes having starting threads;
	Set_start_TID(start_TID, N);
	//cout << "first TID info for each plane is \n";
	//for (int i = 0; i < (2 * N - 1) * 2; i++)
	//cout << start_TID[i] << endl;

	time_used_total = omp_get_wtime();
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
			sweep_ase();
		else if(sweepfun == "eas")
			sweep_eas();
		else if(sweepfun == "esa")
			sweep_esa();
		else if(sweepfun == "sae")
			sweep_sae();
		else if(sweepfun == "sea")
			sweep_sea();

		real sweep_end_time = omp_get_wtime();
		time_used_sweep += sweep_end_time - sweep_begin_time;

		real max = fabs(phi[0] - phi0[0]) / phi0[0];
		for (int i = 1; i < phi_size; i++)
			if (fabs(phi[i] - phi0[i]) / phi0[i] > max)
				max = fabs(phi[i] - phi0[i]) / phi0[i];
		err_phi = max;
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
					//cout << TID << endl;
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
									for(int y0 = blockID_y * block_size; y0 < blockID_y * block_size + block_size; y0++)
									{
										const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
										const int y_local = y_global % block_size;
										const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];

										c = (muDelta * bd_info_x[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, x_local)] +
												etaDelta * bd_info_y[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, y_local)] +
												xiDelta * bd_info_z[x_local * block_size + y_local] + Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e])
																	/ (SigT[m * n_eg + e] + sum);
										phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c;
										bd_info_x[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, x_local)] = 2.0 * c -
												bd_info_x[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, x_local)];
										bd_info_y[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, y_local)] = 2.0 * c -
												bd_info_y[Get_index(N, block_size, blockID_x, blockID_y, blockID_z, z_local, y_local)];
										bd_info_z[x_local * block_size + y_local] = 2.0 * c - bd_info_z[x_local * block_size + y_local];
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

// ****************************************************** division of KBA and without KBA ******************************************************

void Solver::sweep_ase()
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
#pragma omp parallel num_threads(nTs)
	{
		//Here are the thread-private variables.
		real c[n_eg];// psi
		real ch[n_eg]; //horizontal edge flux
		real cv[totNFM_y * n_eg]; //vertical edge flux
		real cz[totNFM_x * totNFM_y * n_eg]; //edge flux in the z direction
		bool forward_x, forward_y, forward_z;

		//Every thread has its own phi_private and initials it to zero.
		real phi_private[phi_size];
		SetValue(phi_private, phi_size, 0.0);

#pragma omp for collapse(2)
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			//Angle loop
			for(int a = 0; a < N_A; a++)
			{
				//------------------------------------------------run in serial------------------------------------------------
				DetermineDir(o, forward_z, forward_x, forward_y);
				const real muDelta = 2.0 * mu[a] / Delta_y;
				const real etaDelta = 2.0 * eta[a] / Delta_x;
				const real xiDelta = 2.0 * xi[a] / Delta_z;
				const real sum = muDelta + etaDelta + xiDelta;

				//spatial loop : z
				SetValue(cz, totNFM_x * totNFM_y * n_eg, 0.0); //vacuum boundary
				for(int z0 = 0; z0 < totNFM_z; z0++)
				{
					const int z = forward_z ? z0 : totNFM_z - 1 - z0;
					SetValue(cv, totNFM_y * n_eg, 0.0);
					//Space loop:x
					for(int x0 = 0; x0 < totNFM_x; x0++)
					{
						const int x = forward_x ? x0 : totNFM_x - 1 - x0;
						SetValue(ch, n_eg, 0.0);
						//Space loop:y
						for(int y0 = 0; y0 < totNFM_y; y0++)
						{
							const int y = forward_y ? y0 : totNFM_y - 1 - y0;
							const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
							//Energy loop
							for(int e = 0; e < n_eg; ++e)
							{
								c[e] = (muDelta * ch[e] + etaDelta * cv[y * n_eg + e] + xiDelta * cz[x * totNFM_y * n_eg + y * n_eg + e] + Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e])
																																																																																																										/ (SigT[m * n_eg + e] + sum);
								phi_private[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c[e];
								ch[e] = 2.0 * c[e] - ch[e];
								cv[y * n_eg + e] = 2.0 * c[e] - cv[y * n_eg + e];
								cz[x * totNFM_y * n_eg + y * n_eg + e] = 2.0 * c[e] - cz[x * totNFM_y * n_eg + y * n_eg + e];
							}
						}
					}
				}
			}
		}

#pragma omp critical
		for(int i = 0; i < phi_size; i++)
			phi[i] += phi_private[i];
	}//End of the parallel region.
}

void Solver::sweep_eas()
{
	const real weight = 0.5 * M_PI / (n_a * n_a);
	SetValue(phi, phi_size, 0.0);
#pragma omp parallel num_threads(nTs)
	{
		//Here are the thread-private variables.
		real c;// psi
		real ch; //horizontal edge flux
		real cv[totNFM_y]; //vertical edge flux
		real cz[totNFM_x * totNFM_y]; //edge flux in the z direction
		bool forward_x, forward_y, forward_z;

		//Every thread has its own phi_private and initials it to zero.
		real phi_private[phi_size];
		SetValue(phi_private, phi_size, 0.0);

#pragma omp for collapse(3)
		//Energy loop
		for(int e = 0; e < n_eg; ++e)
		{
			//Octant loop
			for(int o = 0; o < 8; o++)
			{
				//Angle loop
				for(int a = 0; a < N_A; a++)
				{
					//------------------------------------------------run in serial------------------------------------------------
					DetermineDir(o, forward_z, forward_x, forward_y);
					const real muDelta = 2.0 * mu[a] / Delta_y;
					const real etaDelta = 2.0 * eta[a] / Delta_x;
					const real xiDelta = 2.0 * xi[a] / Delta_z;
					const real sum = muDelta + etaDelta + xiDelta;

					//spatial loop : z
					SetValue(cz, totNFM_x * totNFM_y, 0.0); //vacuum boundary
					for(int z0 = 0; z0 < totNFM_z; z0++)
					{
						const int z = forward_z ? z0 : totNFM_z - 1 - z0;
						SetValue(cv, totNFM_y, 0.0);
						//Space loop:x
						for(int x0 = 0; x0 < totNFM_x; x0++)
						{
							const int x = forward_x ? x0 : totNFM_x - 1 - x0;
							ch = 0.0;
							//Space loop:y
							for(int y0 = 0; y0 < totNFM_y; y0++)
							{
								const int y = forward_y ? y0 : totNFM_y - 1 - y0;
								const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
								c = (muDelta * ch + etaDelta * cv[y] + xiDelta * cz[x * totNFM_y + y] + Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e])
												/ (SigT[m * n_eg + e] + sum);
								phi_private[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c;
								ch = 2.0 * c - ch;
								cv[y] = 2.0 * c - cv[y];
								cz[x * totNFM_y + y] = 2.0 * c - cz[x * totNFM_y + y];
							}
						}
					}
				}
			}
		}

#pragma omp critical
		for(int i = 0; i < phi_size; i++)
			phi[i] += phi_private[i];
	}//End of the parallel region.
}


void Solver::sweep_esa()
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
#pragma omp parallel num_threads(nTs)
	{
		//Here are the thread-private variables.
		real c[N_A];// psi
		real ch[N_A]; //horizontal edge flux
		real cv[totNFM_y * N_A]; //vertical edge flux
		real cz[totNFM_x * totNFM_y * N_A]; //edge flux in the z direction
		bool forward_x, forward_y, forward_z;

		//Every thread has its own phi_private and initials it to zero.
		real phi_private[phi_size];
		SetValue(phi_private, phi_size, 0.0);

#pragma omp for collapse(2)
		//Energy loop
		for(int e = 0; e < n_eg; ++e)
		{
			//Octant loop
			for(int o = 0; o < 8; o++)
			{
				//------------------------------------------------run in serial------------------------------------------------
				DetermineDir(o, forward_z, forward_x, forward_y);
				//spatial loop : z
				SetValue(cz, totNFM_x * totNFM_y * N_A, 0.0); //vacuum boundary
				for(int z0 = 0; z0 < totNFM_z; z0++)
				{
					const int z = forward_z ? z0 : totNFM_z - 1 - z0;
					SetValue(cv, totNFM_y * N_A, 0.0);
					//Space loop:x
					for(int x0 = 0; x0 < totNFM_x; x0++)
					{
						const int x = forward_x ? x0 : totNFM_x - 1 - x0;
						SetValue(ch, N_A, 0.0);
						//Space loop:y
						for(int y0 = 0; y0 < totNFM_y; y0++)
						{
							const int y = forward_y ? y0 : totNFM_y - 1 - y0;
							const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
							//Angle loop
							for(int a = 0; a < N_A; a++)
							{
								const real muDelta = 2.0 * mu[a] / Delta_y;
								const real etaDelta = 2.0 * eta[a] / Delta_x;
								const real xiDelta = 2.0 * xi[a] / Delta_z;
								const real sum = muDelta + etaDelta + xiDelta;
								c[a] = (muDelta * ch[a] + etaDelta * cv[y * N_A + a] + xiDelta * cz[x * totNFM_y * N_A + y * N_A + a] + Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e])
																																																																																																																						/ (SigT[m * n_eg + e] + sum);
								phi_private[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c[a];
								ch[a] = 2.0 * c[a] - ch[a];
								cv[y * N_A + a] = 2.0 * c[a] - cv[y * N_A + a];
								cz[x * totNFM_y * N_A + y * N_A + a] = 2.0 * c[a] - cz[x * totNFM_y * N_A + y * N_A + a];
							}
						}
					}
				}
			}
		}

#pragma omp critical
		for(int i = 0; i < phi_size; i++)
			phi[i] += phi_private[i];
	}//End of the parallel region.
}

void Solver::sweep_sae()
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
#pragma omp parallel num_threads(nTs)
	{
		//Here are the thread-private variables.
		real c[N_A * n_eg];// psi
		real ch[N_A * n_eg]; //horizontal edge flux
		real cv[totNFM_y * N_A * n_eg]; //vertical edge flux
		real cz[totNFM_x * totNFM_y * N_A * n_eg]; //edge flux in the z direction
		bool forward_x, forward_y, forward_z;

		//Every thread has its own phi_private and initials it to zero.
		real phi_private[phi_size];
		SetValue(phi_private, phi_size, 0.0);

#pragma omp for
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			//------------------------------------------------run in serial------------------------------------------------
			DetermineDir(o, forward_z, forward_x, forward_y);
			//spatial loop : z
			SetValue(cz, totNFM_x * totNFM_y * N_A * n_eg, 0.0); //vacuum boundary
			for(int z0 = 0; z0 < totNFM_z; z0++)
			{
				const int z = forward_z ? z0 : totNFM_z - 1 - z0;
				SetValue(cv, totNFM_y * N_A * n_eg, 0.0);
				//Space loop:x
				for(int x0 = 0; x0 < totNFM_x; x0++)
				{
					const int x = forward_x ? x0 : totNFM_x - 1 - x0;
					SetValue(ch, N_A * n_eg, 0.0);
					//Space loop:y
					for(int y0 = 0; y0 < totNFM_y; y0++)
					{
						const int y = forward_y ? y0 : totNFM_y - 1 - y0;
						const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
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
								c[a * n_eg + e] = (muDelta * ch[a * n_eg + e] + etaDelta * cv[y * N_A * n_eg + a * n_eg + e] + xiDelta * cz[x * totNFM_y * N_A * n_eg + y * N_A * n_eg + a * n_eg + e] +
										Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e])	/ (SigT[m * n_eg + e] + sum);
								phi_private[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c[a * n_eg + e];
								ch[a * n_eg + e] = 2.0 * c[a * n_eg + e] - ch[a * n_eg + e];
								cv[y * N_A * n_eg + a * n_eg + e] = 2.0 * c[a * n_eg + e] - cv[y * N_A * n_eg + a * n_eg + e];
								cz[x * totNFM_y * N_A * n_eg + y * N_A * n_eg + a * n_eg + e] = 2.0 * c[a * n_eg + e] - cz[x * totNFM_y * N_A * n_eg + y * N_A * n_eg + a * n_eg + e];
							}
						}
					}
				}
			}
		}

#pragma omp critical
		for(int i = 0; i < phi_size; i++)
			phi[i] += phi_private[i];
	}//End of the parallel region.
}

void Solver::sweep_sea()
{
	const real weight = 0.5 * M_PI / N_A;
	SetValue(phi, phi_size, 0.0);
#pragma omp parallel num_threads(nTs)
	{
		//Here are the thread-private variables.
		real c[n_eg * N_A];// psi
		real ch[n_eg * N_A]; //horizontal edge flux
		real cv[totNFM_y * n_eg * N_A ]; //vertical edge flux
		real cz[totNFM_x * totNFM_y * n_eg * N_A]; //edge flux in the z direction
		bool forward_x, forward_y, forward_z;

		//Every thread has its own phi_private and initials it to zero.
		real phi_private[phi_size];
		SetValue(phi_private, phi_size, 0.0);

#pragma omp for
		//Octant loop
		for(int o = 0; o < 8; o++)
		{
			DetermineDir(o, forward_z, forward_x, forward_y);

			//------------------------------------------------run in serial------------------------------------------------

			//spatial loop : z
			SetValue(cz, totNFM_x * totNFM_y * n_eg * N_A , 0.0); //vacuum boundary
			for(int z0 = 0; z0 < totNFM_z; z0++)
			{
				const int z = forward_z ? z0 : totNFM_z - 1 - z0;
				SetValue(cv, totNFM_y * n_eg * N_A , 0.0);
				//Space loop:x
				for(int x0 = 0; x0 < totNFM_x; x0++)
				{
					const int x = forward_x ? x0 : totNFM_x - 1 - x0;
					SetValue(ch, n_eg * N_A, 0.0);
					//Space loop:y
					for(int y0 = 0; y0 < totNFM_y; y0++)
					{
						const int y = forward_y ? y0 : totNFM_y - 1 - y0;
						const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
						//Energy loop
						for(int e = 0; e < n_eg; ++e)
						{
							//Angle loop
							for(int a = 0; a < N_A; a++)
							{
								const real muDelta = 2.0 * mu[a] / Delta_y;
								const real etaDelta = 2.0 * eta[a] / Delta_x;
								const real xiDelta = 2.0 * xi[a] / Delta_z;
								const real sum = muDelta + etaDelta + xiDelta;
								c[e * N_A + a] = (muDelta * ch[e * N_A + a] + etaDelta * cv[y * n_eg * N_A + e * N_A + a] + xiDelta * cz[x * totNFM_y * n_eg * N_A + y * n_eg * N_A + e * N_A + a] +
										Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e])	/ (SigT[m * n_eg + e] + sum);
								phi_private[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c[e * N_A + a];
								ch[e * N_A + a] = 2.0 * c[e * N_A + a] - ch[e * N_A + a];
								cv[y * n_eg * N_A + e * N_A + a] = 2.0 * c[e * N_A + a] - cv[y * n_eg * N_A + e * N_A + a];
								cz[x * totNFM_y * n_eg * N_A + y * n_eg * N_A + e * N_A + a] = 2.0 * c[e * N_A + a] - cz[x * totNFM_y * n_eg * N_A + y * n_eg * N_A + e * N_A + a];
							}
						}
					}
				}
			}
		}

#pragma omp critical
		for(int i = 0; i < phi_size; i++)
			phi[i] += phi_private[i];
	}//End of the parallel region.
}

real Solver::get_sweeptime(){
	return time_used_sweep;
}

real Solver::get_totaltime(){
	return time_used_total;
}
