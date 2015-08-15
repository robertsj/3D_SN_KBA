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
	Solver test(4, 5, 4, 6, 2, 5);
	test.Calculate("aes", 4);
	test.Calculate("ase", 4);
	test.Calculate("eas", 4);
	test.Calculate("esa", 4);
	test.Calculate("sae", 4);
	test.Calculate("sea", 4);
	return 0;
}



