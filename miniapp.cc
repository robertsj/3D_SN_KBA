#include <iostream>
#include <cmath>
#include <time.h>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <cassert>
#include <string>

#include "miniapp.hh"
#include "auxiliary_function.hh"

using namespace std;

typedef std::vector<real> v_dbl;
#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

//----------------------------------------------------------------------------//
Solver::Solver(std::string order_in, int n_eg_in, int n_a_in, int cm_xy_in, int fm_xy_in,
    int cm_z_in, int fm_z_in, int upscatter_in, int iter_in,
    int xbs, int ybs, int zbs, int nt)
  : order(order_in)
  , n_eg(n_eg_in)
  , n_a(n_a_in)
  , cm_xy(cm_xy_in)
  , fm_xy(fm_xy_in)
  , cm_z(cm_z_in)
  , fm_z(fm_z_in)
  , upscatter(upscatter_in)
  , iter(iter_in)
  , blocksize_xy(xbs)
  , blocksize_z(zbs)
  , nTs(nt)
  , mesh(fm_xy, cm_xy, fm_z, cm_z, xbs, xbs, zbs)
{

  cout << "input parameters \n" << "# of energy groups: " << n_eg << "\n"
      << "# of angles in each octant: " << n_a << "\n"
      << "# of coarse cells along x & y axes: " << cm_xy << " \n"
      << "# of fine cells in one coarse cell along x & y axes: "<< fm_xy << "\n"
      << "# of coarse cells along z axis: " << cm_z << " \n"
      << "# of fine cells in one coarse cell along z axis: "<< fm_z << "\n"
      << "# of upscatter energy groups: " << upscatter << "\n"
      << "# of running iterations: " << iter << "\n"
      << "xy block size: " << xbs << endl
      << " z block size: " << zbs << endl;

  totNFM_x = cm_xy * fm_xy;
  totNFM_y = cm_xy * fm_xy;
  totNFM_z = cm_z * fm_z;
  phi_size = totNFM_z * totNFM_x * totNFM_y  * n_eg;

  phi = new real [phi_size];
  Q = new real [phi_size];

  //Since the core mesh and the fine mesh is uniformly distributed over the total area, so we have the following Delta_x,Delta_y.
  //In the case the total length in X and Y direction is 1.0.
  Delta_x = 1.0 / totNFM_x;
  Delta_y = 1.0 / totNFM_y;
  Delta_z = 1.0 / totNFM_z;

  n_m = cm_xy * cm_xy * cm_z;
  //Let each part of the core mesh contain a random unique material
  srand(time(0));
  for(int i = 0; i < n_m; i++)
    RegMat.push_back(i);
  random_shuffle(RegMat.begin(), RegMat.end());

  //allocate and initiate fmmid, fine mesh (cell) material ID
  fmmid = new int [totNFM_x * totNFM_y * totNFM_z];
  for (int z = 0; z < cm_z; z ++)
    for(int x = 0; x < cm_xy; x ++)
      for (int y = 0; y < cm_xy; y ++)
      {
        int m = RegMat[z * cm_xy * cm_xy + x * cm_xy + y];
        for (int zfm = z * fm_z; zfm < (z + 1) * fm_z; zfm++)
          for(int xfm = x * fm_xy; xfm < (x + 1) * fm_xy; xfm++)
            for(int yfm = y * fm_xy; yfm < (y + 1) * fm_xy; yfm++)
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

//----------------------------------------------------------------------------//
void Solver::get_quadrature()
{
  /*  This provides an easy but sort of realistic quadrature defined only
   *  by the numer of angles per octant.
   */
  mu.resize(n_a);
  eta.resize(n_a);
  xi.resize(n_a);
  wt.resize(n_a);
  for (int i = 0; i < n_a; ++i)
  {
    mu[i] = 1.0/(i+0.5/i/n_a);
    int j = n_a - i - 1;
    eta[i] = 1.0/(j+0.5/j/n_a);
    xi[i] = std::sqrt(1.0-mu[i]*mu[i]-eta[i]*eta[i]);
    wt[i] = 1.5707963267948966 / n_a;
  }
}

void Solver::Calculate()
{
  time_used_total = omp_get_wtime();
  assert(order == "aes" || order == "ase" || order == "eas" ||
         order == "esa" || order == "sae" || order == "sea" ||
         order == "eas_mod");
  //nTs = nTs_in;
  N = sqrt(nTs);
  assert(totNFM_x % N == 0);
  assert(blocksize_z <= totNFM_z);
  cout << "# of threads is " << nTs << endl;
  blocksize_xy = totNFM_x / N;
  //blocksize_z = blocksize_z_in;
  n_b = totNFM_z / blocksize_z;
  n_p = 2 * N - 2 + n_b;
  remain = totNFM_z % blocksize_z;

  int start_TID[(2 * N - 1) * 2];// start_TID -- TID info for planes having starting threads;
  Set_start_TID(start_TID, N);

  //Initial guess of phi, 0.0
  SetValue(phi, phi_size, 0.0);

  time_used_sweep = 0.0;

  //convergence parameters
  real eps_phi = 1.0e-5;
  real err_phi = 1.0;
  int it = 0;

  while(it < iter)
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

    if(order == "aes")
      sweep_aes(start_TID);
    else if(order == "ase")
      sweep_ase(start_TID);
    else if(order == "eas")
      sweep_eas(start_TID);
    else if(order == "esa")
      sweep_esa(start_TID);
    else if(order == "sae")
      sweep_sae(start_TID);
    else if(order == "sea")
      sweep_sea(start_TID);
    else if(order == "eas_mod")
      sweep_eas_mod(start_TID);

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
    cout << "time used in " << order << " sweep of " << it << " iterations is " << time_used_sweep << " s" << endl;
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
    cout << "\ndo not converge and time used in "<< order << " sweep in " << iter << " iterations is " << time_used_sweep << " s" << endl;
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

real Solver::get_sweeptime(){
  return time_used_sweep;
}

real Solver::get_totaltime(){
  return time_used_total;
}
