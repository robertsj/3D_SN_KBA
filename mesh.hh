#ifndef MESH_HH_
#define MESH_HH_

#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

#include <vector>

/**
 *  Class to automate the construction of a mesh, the random assignment of
 *  material indices, and the decomposition of the mesh for cache-optimal
 *  threading.
 */
class Mesh
{

public:

  typedef std::vector<real> vec_dbl;
  typedef std::vector<int> vec_int;
  typedef std::vector<vec_int> vec2_int;
  typedef std::vector<vec2_int> vec3_int;
  typedef std::vector<vec3_int> vec4_int;

  enum AXIS
  {
    X, Y, Z
  };

  Mesh(int fm_xy, int cm_xy, int fm_z, int cm_z, int bs_x, int bs_y, int bs_z);
  void print_blocks();

  // go from (i, j, k) to cardinal index
  int ijk_to_cell(int i, int j, int k) const
  {
    return nx*ny*k + nx*j + i;
  }

  // for a given octant, are we going forward in x, y, and z?
  void forward(const int o, int &xflag, int &yflag, int &zflag) const
  {
    xflag = o % 4 == 0 or o % 4 == 3;
    yflag = o % 4 < 2;
    zflag = o < 4;
  }

  // given local i, j, or k, get the actual i, j, and k.  the point is that
  // each block can be looped over using, e.g., "for i in range(0, xbs[i_b])".
  int get_block_i(const int o, int i_b, int i) const
  {
    return (o % 4 == 0 or o % 4 == 3) ? i+i_s[i_b] : i_s[i_b]+xbs[i_b]-i-1;
  }
  int get_block_j(const int o, int j_b, int j) const
  {
    return (o % 4 < 2) ? j+j_s[j_b] : j_s[j_b]+ybs[j_b]-j-1;
  }
  int get_block_k(const int o, int k_b, int k) const
  {
    return (o < 4) ? k+k_s[k_b] : k_s[k_b]+zbs[k_b]-k-1;
  }

  int nx, ny, nz, ncell;            // # of cells along x, y, and z and total #
  vec_dbl dx, dy, dz;               // cell widths
  int nmat;                         // number of materials
  vec_int mat_id;                   // material id for each cell
  int nbx, nby, nbz;                // number of blocks along x, y, and z
  vec_int xbs, ybs, zbs;            // block sizes
  vec_int i_s, j_s, k_s;            // starting cell indices for each block
  vec4_int planes;                  // planes with indices
};

#endif
