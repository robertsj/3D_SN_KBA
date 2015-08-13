/*
 * main.cc
 *
 *  Created on: Aug 13, 2015
 *      Author: kevin
 */
#include <iostream>
#include "miniapp.hh"
#include <string>

int main()
{
	Solver test(10, 1, 2, 2, 0, 100);
	test.Calculate("aes", 4);
	return 0;
}



