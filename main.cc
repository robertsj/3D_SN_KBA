/*
 * main.cc
 *
 *  Created on: Aug 13, 2015
 *      Author: kevin
 */
#include "miniapp.hh"
#include <string>
using namespace std;
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>


int main()
{
  //          ng  a  cm  fm cmz fmz up iter
	Solver test(1, 8, 1, 240, 1, 240, 0, 5);
	
	int nt[] = {1, 4, 9, 16};
    for (int i = 0; i < 4; ++i)
	    test.Calculate("esa",  nt[i], 1);
  //test.Calculate("sea", 16, 1);

	return 0;
}



