/*
 * auxiliary_function.hh
 *
 *  Created on: Aug 12, 2015
 *      Author: kevin
 */

#ifndef AUXILIARY_FUNCTION_HH_
#define AUXILIARY_FUNCTION_HH_

#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

void SetValue(real *array, int length, real value);

void DetermineDir(int o, bool &forward_z, bool &forward_x, bool &forward_y);


void copy(real *left, int lstart, real *right, int rstart, int length);


void Set_start_TID(int start_TID[], int N);


void Find_start_plane(int start_TID[], int N, int TID, int &start_plane);


void Set_block_ID_xy(int start_TID[], int N, int start_plane, int TID, int &block_x, int &block_y);

//inline int Get_index(int N, int block_size, int blockID_x, int blockID_y, int blockID_z, int z_local, int x_or_y_local);
#endif /* AUXILIARY_FUNCTION_HH_ */
