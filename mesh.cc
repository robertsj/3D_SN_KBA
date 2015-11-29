#include "mesh.hh"
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <algorithm>
#include "auxiliary_function.hh"

//----------------------------------------------------------------------------//
Mesh::Mesh(int fm_xy, int cm_xy,
           int fm_z, int cm_z,
           int bs_x, int bs_y, int bs_z)
  : nx(fm_xy * cm_xy)
  , ny(fm_xy * cm_xy)
  , nz(fm_z  * cm_z)
  , ncell(nx*ny*nz)
  , dx(nx, 10.0/nx)
  , dy(ny, 10.0/ny)
  , dz(nz, 10.0/nz)
  , nmat(cm_xy*cm_xy*cm_z)
  , mat_id(ncell, 0)
  , nbx(nx / bs_x + (1 ? nx % bs_x : 0))
  , nby(ny / bs_y + (1 ? ny % bs_y : 0))
  , nbz(nz / bs_z + (1 ? nz % bs_z : 0))
  , xbs(nbx, bs_x)
  , ybs(nby, bs_y)
  , zbs(nbz, bs_z)
  , i_s(nbx)
  , j_s(nby)
  , k_s(nbz)
  , planes(8)
{
  printf("mesh details: \n");
  printf("    nx = %i\n", nx);
  printf("    ny = %i\n", ny);
  printf("    nz = %i\n", nz);
  printf("  nmat = %i\n", nmat);


  // fix size of last blocks block sizes
  xbs[nbx-1] = nx - bs_x * (nbx - 1);
  ybs[nby-1] = ny - bs_y * (nby - 1);
  zbs[nbz-1] = nz - bs_z * (nbz - 1);

  // randomize coarse mesh materials
  std::srand(time(0));
  vec_int reg_mat(nmat);
  for(int i = 0; i < nmat; ++i)
    reg_mat[i] = i;
  std::random_shuffle(mat_id.begin(), mat_id.end());
  for (int z = 0; z < cm_z; ++z)
  {
    for(int y = 0; y < cm_xy; ++y)
    {
      for (int x = 0; x < cm_xy; ++x)
      {
        int m = reg_mat[z * cm_xy * cm_xy + y * cm_xy + x];
        // assign this material to all cells in this coarse cell
        for (int k = z * fm_z; k < (z + 1) * fm_z; ++k)
          for(int j = y * fm_xy; j < (y + 1) * fm_xy; ++j)
            for(int i = x * fm_xy; i < (x + 1) * fm_xy; ++i)
              mat_id[ijk_to_cell(i, j, k)] = m;
      }
    }
  }

  // define planes
  int np = nbx + nby + nbz - 2;
  for (int o = 0; o < 8; ++o)
  {
    planes[o].resize(np);
    for (int p = 0; p < np; ++p)
    {
      for (int ii = 0; ii < nbx; ++ii)
      {
        int i = (o % 4 == 0 or o % 4 == 3) ? ii : nbx - ii - 1;
        for (int jj = 0; jj < nby; ++jj)
        {
          int j = (o % 4 < 2) ? jj : nby - jj - 1;
          for (int kk = 0; kk < nbz; ++kk)
          {
            int k =  (o < 4) ? kk : nbz - kk - 1;
            if (ii + jj + kk == p)
            {
              vec_int ijk(3);
              ijk[0] = i;
              ijk[1] = j;
              ijk[2] = k;
              planes[o][p].push_back(ijk);
            }
          }
        }
      }
    }
  } // end octant

  // block indices
  i_s.resize(nbx, 0);
  j_s.resize(nby, 0);
  k_s.resize(nbz, 0);

  for (int i = 0; i < nbx; ++i)
    for (int ii = 0; ii < i; ++ii)
      i_s[i] += xbs[ii];

  for (int j = 0; j < nby; ++j)
    for (int jj = 0; jj < j; ++jj)
      j_s[j] += ybs[jj];

  for (int k = 0; k < nbz; ++k)
    for (int kk = 0; kk < k; ++kk)
      k_s[k] += zbs[kk];
}



//----------------------------------------------------------------------------//
void Mesh::print_blocks()
{
  print_dline();
  for (int o = 0; o < 8; ++o)
  {
    print_dline();
    for (int p = 0; p < planes[o].size(); ++p)
    {
      vec2_int plane = planes[o][p];
      for (int b = 0; b < plane.size(); ++b)
      {
        vec_int block = plane[b];
        int i0=get_block_i(o, block[X], 0);
        int i1=get_block_i(o, block[X], xbs[block[X]]-1);
        int j0=get_block_j(o, block[Y], 0);
        int j1=get_block_j(o, block[Y], ybs[block[Y]]-1);
        int k0=get_block_k(o, block[Z], 0);
        int k1=get_block_k(o, block[Z], zbs[block[Z]]-1);
        printf("octant=%i, plane=%i, block=%i, i=%i-->%i, j=%i-->%i, k=%i-->%i\n",
            o, p, b, i0, i1, j0, j1, k0, k1);
      }
      print_line();
    }
    print_dline();
  }

  for (int i = 0; i < k_s.size(); ++i)
  {
    printf("%i %i\n", i, zbs[i]);
  }
}
