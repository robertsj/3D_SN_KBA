#include "miniapp.hh"
#include "auxiliary_function.hh"

void Solver::sweep_eas_mod(int start_TID[])
{
  SetValue(phi, phi_size, 0.0);

  // use, e.g., ybs[0] because only ybs[-1] could be smaller
  // x-directed  Z         Y         z size       y size
  real bd_info_x[mesh.nbz][mesh.nby][mesh.zbs[0]][mesh.ybs[0]];
  // y-directed  Z         X         z size       x size
  real bd_info_y[mesh.nbz][mesh.nbx][mesh.zbs[0]][mesh.xbs[0]];
  // z-directed  Y         X         y size       x size
  real bd_info_z[mesh.nby][mesh.nbx][mesh.ybs[0]][mesh.xbs[0]];

  real muDelta[n_a], etaDelta[n_a], xiDelta[n_a], sum[n_a];
  for(int a = 0; a < n_a; a++)
  {
    muDelta[a]  = 2.0 * mu[a] / Delta_y;
    etaDelta[a] = 2.0 * eta[a] / Delta_x;
    xiDelta[a]  = 2.0 * xi[a] / Delta_z;
    sum[a]      = muDelta[a] + etaDelta[a] + xiDelta[a];
  }

  // material index
  const int *mat_id = &mesh.mat_id[0];

  // GROUPS
  for (int g = 0; g < n_eg; ++g)
  {
    // OCTANTS
    for (int o = 0; o < 8; o++)
    {
      // ANGLES
      for(int a = 0; a < n_a; a++)
      {
        // set all cell edge fluxes to zero.  this is consistent with vacuum.
        SetValue(&bd_info_x[0][0][0][0], mesh.nbz*mesh.nby*mesh.zbs[0]*mesh.ybs[0], 0.0);
        SetValue(&bd_info_y[0][0][0][0], mesh.nbz*mesh.nbx*mesh.zbs[0]*mesh.xbs[0], 0.0);
        SetValue(&bd_info_z[0][0][0][0], mesh.nby*mesh.nbx*mesh.ybs[0]*mesh.xbs[0], 0.0);

#pragma omp parallel num_threads(nTs)
        {
          int TID = omp_get_thread_num();
          // PLANES
          for (int p = 0; p < mesh.planes[o].size(); ++p)
          {
            // BLOCKS
            for (int b = 0; b < mesh.planes[o][p].size(); ++b)
            {
              // SKIP THIS BLOCK IF NOT MINE
              if (b % nTs != TID)
                continue;
              // block id's
              int i_b = mesh.planes[o][p][b][mesh.X];
              int j_b = mesh.planes[o][p][b][mesh.Y];
              int k_b = mesh.planes[o][p][b][mesh.Z];
              // Z
              for (int kk = 0; kk < mesh.zbs[k_b]; ++kk)
              {
                int k = mesh.get_block_k(o, k_b, kk);
                // Y
                for (int jj = 0; jj < mesh.ybs[j_b]; ++jj)
                {
                  int j = mesh.get_block_j(o, j_b, jj);
                  // X
                  for (int ii = 0; ii < mesh.xbs[i_b]; ++ii)
                  {
                    int i = mesh.get_block_i(o, i_b, ii);
                    int cell = mesh.ijk_to_cell(i, j, k);
                    int m = mat_id[cell];
                    double qq=1.0, sigt=1.0;
                    //========================================================//
                    // THE WORK!
                    {
                      real c = ( muDelta[a] * bd_info_x[k_b][j_b][kk][jj] +
                                etaDelta[a] * bd_info_y[k_b][i_b][kk][ii] +
                                 xiDelta[a] * bd_info_z[j_b][i_b][jj][ii] +
                                 qq) / (sigt + sum[a]);
                      phi[cell] += wt[a] * c;
                      bd_info_x[k_b][j_b][kk][jj] = 2.0*c - bd_info_x[k_b][j_b][kk][jj];
                      bd_info_y[k_b][i_b][kk][ii] = 2.0*c - bd_info_y[k_b][i_b][kk][ii];
                      bd_info_z[j_b][i_b][jj][ii] = 2.0*c - bd_info_z[j_b][i_b][jj][ii];
                    }
                    //========================================================//

                  } // I
                } // J
              } // K
            } // BLOCKS
#pragma omp barrier
          } // PLANES
        } // PARALLEL
      } // ANGLES
    } // OCTANTS
  } // GROUPS
}
