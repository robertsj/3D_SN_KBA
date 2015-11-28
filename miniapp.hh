/*
 * miniapp.hh
 *
 *  Created on: Nov 4, 2014
 *      Author: kevin
 */
#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

#ifndef MINIAPP_HH_
#define MINIAPP_HH_
#include <string>
#include <iostream>
#include <cmath>
#include <time.h>
#include <omp.h>

using namespace std;
#include <vector>
typedef std::vector<real> v_dbl;

//Assume X,Y directions are equally divided by the cm core mesh,each core mesh is equally divided by the fm fine mesh.
//Assume that the length in X,Y direction are 1.0.
//The cross section data is set inside the class.
//Assume there is one Source=1.0 in each core mesh.

class Solver
{

private:

  /* n_eg -- # of energy groups;
   * n_a -- in each octant, there are n_a * n_a = N_A angles;
   * n_m -- # of materials;
   * cm -- # of coarse cells in one direction;
   * fm -- # of fine cells in a coarse cell along a direction;
   * for 3D, the spatial size is (cm * fm) ^ 3, i.e., a cubic geometry;
   * upscatter -- # of energy groups consider upscarttering;
   * iter -- # of running iteration;
   * totNFM_x, _y and _z -- total number of fine cells in each direction direction, which equals to cm * fm;
   * nTs -- # of threads used, in 3D KBA, must be i ^ 2, e.g., 1, 4, 9, 16...;
   * N -- sqrt of nTs, also means # of blocks in x & y directions;
   * block_size -- block size along each direction, = totNFM_x or _y or _z / N;
   * sweepfun -- used to choose sweep function, aes, ase, eas, esa, sae, sea;
   * mu, eta, xi -- cosines of neutron direction with x, y, z axes, respectively;
   * phi -- scaler flux, [Z][X][Y][E];
   * phi_size -- = totNFM_z * totNFM_x * totNFM_y * n_eg;
   * Delta_x, _y and _z -- width of a fine cell in each direction;
   * RegMat[k][i][j] -- material ID in coarse cell (i, j, k), material ID is randomly assigned;
   * note the ordinates sequence (ZXY) to realize unit strid, which is the same for all the space-dependent variables, e.g., phi, Q;
   * SigT[i][j] -- total cross section of material i in energy group j;
   * SigS[i][j][k] -- scattering cross section of material i,scattering from energy group k to j;
   * If k > j, SigS = 0, except the energy group consider the upscattering;
   * begin_time:used to calculate time used;
   * Q -- emission density, has phase space same as phi, i.e., [Z][X][Y][E];
   * fmmid[z][x][y] -- material ID in fine cell [z][x][y];
   * n_b -- # of blocks in z direction, = totNFM_z / blocksize_z
   * n_p -- # of diagonal planes in block sweep
   * remain -- remain fine cells along z axis after block sweep, = totNFM_z % blocksize_z

   * */

  int n_eg, n_a, N_A, n_m, cm_xy, fm_xy, cm_z, fm_z, phi_size, upscatter, iter, totNFM_x, totNFM_y, totNFM_z, nTs, N, blocksize_xy;
  int blocksize_z, n_b, n_p, remain;
  v_dbl mu, eta, xi;
  real* phi;
  real Delta_x, Delta_y, Delta_z;
  vector<int> RegMat;
  real* SigT;
  real*** SigS;
  real* Q;
  int* fmmid;
  real time_used_sweep;
  real time_used_total;

public:
  Solver(int n_eg_in, int n_a_in, int cm_xy_in, int fm_xy_in, int cm_z_in, int fm_z_in, int upscatter_in, int iter_in);
  ~Solver();
  void Calculate(string sweepfun, int nTs_in, int blocksize_z_in);
  void get_quadrature();
  void sweep_aes(int start_TID[]);
  void sweep_ase(int start_TID[]);
  void sweep_eas(int start_TID[]);
  void sweep_esa(int start_TID[]);
  void sweep_sae(int start_TID[]);
  void sweep_sea(int start_TID[]);
  real get_sweeptime();
  real get_totaltime();
};


inline int Get_index(int N, int n_b, int blocksize_z, int blocksize_xy, int sweep,
    int blockID_x, int blockID_y, int blockID_z, int z_local, int x_or_y_local)
{
  return blockID_x * (N + 1) * n_b * blocksize_z * blocksize_xy * sweep + blockID_y * n_b * blocksize_z * blocksize_xy * sweep +
      blockID_z * blocksize_z * blocksize_xy * sweep + z_local * blocksize_xy * sweep + x_or_y_local * sweep;
}


#endif /* MINIAPP_HH_ */
