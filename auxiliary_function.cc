/*
 * auxiliary_function.cc
 *
 *  Created on: Aug 12, 2015
 *      Author: kevin
 */
#include "auxiliary_function.hh"
#include <iostream>

#ifdef FLOAT
typedef float real;
#else
typedef double real;
#endif

void SetValue(real *array, int length, real value){
	for (int i = 0; i < length; i++)
		array[i] = value;
}

void DetermineDir(int o, bool &forward_z, bool &forward_x, bool &forward_y)
{
	if (o < 4)
		forward_z = true;
	else
		forward_z = false;
	if(o == 0 || o == 1 || o == 4 || o == 5)
		forward_x = true;
	else
		forward_x = false;
	if(o == 0 || o == 2 || o == 4 || o == 6)
		forward_y = true;
	else
		forward_y = false;
}

void copy(real *left, int lstart, real *right, int rstart, int length)
{
	for (int i = 0; i < length; i++)
		left[lstart + i] = right[rstart + i];
}

void Set_start_TID(int start_TID[], int N)
{
	for (int p = 0; p < N; p++)
	{
		start_TID[p * 2] = p * (p + 1) / 2;
		start_TID[p * 2 + 1] = p * (p + 1) / 2 + p + 1;
	}

	for (int p = N; p < 2 * N - 1; p++)
	{
		start_TID[p * 2] = N * (N - 1) + p - (N - 1) - (2 * N - 1 - p) * (2 * N - 2 - p) / 2;
		start_TID[p * 2 + 1] = start_TID[p * 2] + 2 * N - 1 - p;
	}
}

void Find_start_plane(int start_TID[], int N, int TID, int &start_plane)
{
	for (int p = 0; p < 2 * N - 1; p++)
		if (start_TID[p * 2] <= TID && TID < start_TID[p * 2 + 1])
			start_plane = p;
}

void Set_block_ID_xy(int start_TID[], int N, int start_plane, int TID, int &block_x, int &block_y)
{
	if (start_plane < N)
	{
		block_y = start_plane - (TID - start_TID[start_plane * 2]);
		block_x = TID - start_TID[start_plane * 2];
	}
	else
	{
		block_y = N - 1 - (TID - start_TID[start_plane * 2]);
		block_x = start_plane - (N - 1) + TID - start_TID[start_plane * 2];
	}
}
