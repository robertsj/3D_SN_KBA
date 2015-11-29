#include "miniapp.hh"
#include "auxiliary_function.hh"

void Solver::sweep_ase(int start_TID[])
{
  const real weight = 0.5 * M_PI / n_a;
  SetValue(phi, phi_size, 0.0);
  bool forward_x, forward_y, forward_z;
  //The arrangement of bd_info_x or _y is [blockID_X][blockID_Y][blockID_Z][z][x or y] to realize unit strid
  //one row is wasted along each direction
  real bd_info_x[(N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy * n_eg], bd_info_y[(N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy * n_eg];
  real ch[n_eg], cv[totNFM_y * n_eg], cz[totNFM_x * totNFM_y * n_eg];//used in the remain sweep after block sweep

  real muDelta[n_a], etaDelta[n_a], xiDelta[n_a], sum[n_a];
  for(int a = 0; a < n_a; a++)
  {
    muDelta[a] = 2.0 * mu[a] / Delta_y;
    etaDelta[a] = 2.0 * eta[a] / Delta_x;
    xiDelta[a] = 2.0 * xi[a] / Delta_z;
    sum[a] = muDelta[a] + etaDelta[a] + xiDelta[a];
  }

  //Octant loop
  for(int o = 0; o < 8; o++)
  {
    DetermineDir(o, forward_z, forward_x, forward_y);
    //Angle loop
    for(int a = 0; a < n_a; a++)
    {
      //for simplicity, set all to 0, in fact only the start bd_info need to set 0
      SetValue(bd_info_x, (N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy * n_eg, 0.0);
      SetValue(bd_info_y, (N + 1) * (N + 1) * n_b * blocksize_z * blocksize_xy * n_eg, 0.0);

#pragma omp parallel num_threads(nTs)
      {
        //thread-private variables
        int TID = omp_get_thread_num();
        int start_plane, end_plane;
        Find_start_plane(start_TID, N, TID, start_plane);
        end_plane = start_plane + n_b; //working plane for TID is [start_plane, end_plane);
        //block ID
        int blockID_x, blockID_y, blockID_z;
        Set_block_ID_xy(start_TID, N, start_plane, TID, blockID_x, blockID_y);
        real bd_info_z[blocksize_xy * blocksize_xy * n_eg]; //[X][Y][E]
        SetValue(bd_info_z, blocksize_xy * blocksize_xy * n_eg, 0.0);
        real c[n_eg];

        for(int p = 0; p < n_p; p++)
        {
          if (p >= start_plane && p < end_plane)//select working threads
          {
            blockID_z = p - start_plane;
            //block sweep
            for(int z0 = blockID_z * blocksize_z; z0 < blockID_z * blocksize_z + blocksize_z; z0++) //always use global ordinates
            {
              const int z_global = forward_z ? z0 : totNFM_z - 1 - z0;
              const int z_local = z_global % blocksize_z;
              for(int x0 = blockID_x * blocksize_xy; x0 < blockID_x * blocksize_xy + blocksize_xy; x0++)
              {
                const int x_global = forward_x ? x0 : totNFM_x - 1 - x0;
                const int x_local = x_global % blocksize_xy;
                const int index_x = Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x, blockID_y, blockID_z, z_local, x_local);
                for(int y0 = blockID_y * blocksize_xy; y0 < blockID_y * blocksize_xy + blocksize_xy; y0++)
                {
                  const int y_global = forward_y ? y0 : totNFM_y - 1 - y0;
                  const int y_local = y_global % blocksize_xy;
                  const int m = fmmid[z_global * totNFM_x * totNFM_y + x_global * totNFM_y + y_global];
                  const int index_y = Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x, blockID_y, blockID_z, z_local, y_local);
                  const int index_z = x_local * blocksize_xy * n_eg + y_local * n_eg;
                  //Energy loop
#pragma omp simd
                  for(int e = 0; e < n_eg; ++e)
                  {
                    c[e] = (muDelta[a] * bd_info_x[index_x + e] + etaDelta[a] * bd_info_y[index_y + e] + xiDelta[a] * bd_info_z[index_z + e] +
                        Q[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e]) / (SigT[m * n_eg + e] + sum[a]);
                    phi[z_global * totNFM_x * totNFM_y * n_eg + x_global * totNFM_y * n_eg + y_global * n_eg + e] += weight * c[e];
                    bd_info_x[index_x + e] = 2.0 * c[e] - bd_info_x[index_x + e];
                    bd_info_y[index_y + e] = 2.0 * c[e] - bd_info_y[index_y + e];
                    bd_info_z[index_z + e] = 2.0 * c[e] - bd_info_z[index_z + e];
                  }
                }
              }
            }
            //after block sweep, copy the bd_info_x and _y from the memory of TID into the corresponding memory of (TID + 1)
            copy(bd_info_x, Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x, blockID_y + 1, blockID_z, 0, 0),
                bd_info_x, Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x, blockID_y, blockID_z, 0, 0), blocksize_z * blocksize_xy * n_eg);
            copy(bd_info_y, Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x + 1, blockID_y, blockID_z, 0, 0),
                bd_info_y, Get_index(N, n_b, blocksize_z, blocksize_xy, n_eg, blockID_x, blockID_y, blockID_z, 0, 0), blocksize_z * blocksize_xy * n_eg);
          }
#pragma omp barrier
        }
        if(remain != 0)
        {
          if (forward_x == true && forward_y == true)
            for(int x0 = 0; x0 < blocksize_xy; x0++)
              for(int y0 = 0; y0 < blocksize_xy; y0++)
                for(int e = 0; e < n_eg; e++)
                  cz[(blockID_x * blocksize_xy + x0) * totNFM_y * n_eg + blockID_y * blocksize_xy * n_eg + y0 * n_eg + e] =
                      bd_info_z[x0 * blocksize_xy * n_eg + y0 * n_eg + e];
          else if(forward_x == true && forward_y == false)
            for(int x0 = 0; x0 < blocksize_xy; x0++)
              for(int y0 = 0; y0 < blocksize_xy; y0++)
                for(int e = 0; e < n_eg; e++)
                  cz[(blockID_x * blocksize_xy + x0) * totNFM_y * n_eg + (N - 1 - blockID_y) * blocksize_xy * n_eg + y0 * n_eg + e] =
                      bd_info_z[x0 * blocksize_xy * n_eg + y0 * n_eg + e];
          else if(forward_x == false && forward_y == true)
            for(int x = 0; x < blocksize_xy; x++)
              for(int y = 0; y < blocksize_xy; y++)
                for(int e = 0; e < n_eg; e++)
                  cz[((N - 1 - blockID_x) * blocksize_xy + x) * totNFM_y * n_eg + blockID_y * blocksize_xy * n_eg + y * n_eg + e] =
                      bd_info_z[x * blocksize_xy * n_eg + y * n_eg + e];
          else
            for(int x = 0; x < blocksize_xy; x++)
              for(int y = 0; y < blocksize_xy; y++)
                for(int e = 0; e < n_eg; e++)
                  cz[((N - 1 - blockID_x) * blocksize_xy + x) * totNFM_y * n_eg + (N - 1 - blockID_y) * blocksize_xy * n_eg + y * n_eg + e] =
                      bd_info_z[x * blocksize_xy * n_eg + y * n_eg + e];
        }
      }
      //remain sweep
      if (remain != 0)
      {
        real c[n_eg];
        for(int z0 = 0; z0 < remain; z0++)
        {
          const int z = forward_z ? n_b * blocksize_z + z0 : remain - 1 - z0;
          SetValue(cv, totNFM_y * n_eg, 0.0);
          for(int x0 = 0; x0 < totNFM_x; x0++)
          {
            const int x = forward_x ? x0 : totNFM_x - 1 - x0;
            SetValue(ch, n_eg, 0.0);
            //Space loop:y
            for(int y0 = 0; y0 < totNFM_y; y0++)
            {
              const int y = forward_y ? y0 : totNFM_y - 1 - y0;
              const int m = fmmid[z * totNFM_x * totNFM_y + x * totNFM_y + y];
#pragma omp simd
              for(int e = 0; e < n_eg; e++)
              {
                c[e] = (muDelta[a] * ch[e] + etaDelta[a] * cv[y * n_eg + e] + xiDelta[a] * cz[x * totNFM_y * n_eg + y * n_eg + e] +
                    Q[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e]) / (SigT[m * n_eg + e] + sum[a]);
                phi[z * totNFM_x * totNFM_y * n_eg + x * totNFM_y * n_eg + y * n_eg + e] += weight * c[e];
                ch[e] = 2.0 * c[e] - ch[e];
                cv[y * n_eg + e] = 2.0 * c[e] - cv[y * n_eg + e];
                cz[x * totNFM_y * n_eg + y * n_eg + e] = 2.0 * c[e] - cz[x * totNFM_y * n_eg + y * n_eg + e];
              }
            }
          }
        }
      }
    }
  }
}

