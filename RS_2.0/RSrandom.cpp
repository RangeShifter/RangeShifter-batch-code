/*----------------------------------------------------------------------------
 *	
 *	Copyright (C) 2020 Greta Bocedi, Stephen C.F. Palmer, Justin M.J. Travis, Anne-Kathleen Malchow, Damaris Zurell 
 *	
 *	This file is part of RangeShifter.
 *	
 *	RangeShifter is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *	
 *	RangeShifter is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *	GNU General Public License for more details.
 *	
 *	You should have received a copy of the GNU General Public License
 *	along with RangeShifter. If not, see <https://www.gnu.org/licenses/>.
 *	
 --------------------------------------------------------------------------*/
 
 

#pragma hdrstop

#include "RSrandom.h"
//---------------------------------------------------------------------------


#pragma package(smart_init)

#if RSDEBUG
#include "Parameters.h"
extern paramSim *paramsSim;
#include <fstream>
//ofstream RSRANDOMLOG;
#endif

#if !LINUX_CLUSTER
// set up Mersenne Twister random generator
#if RSWIN64
#if RSDEBUG
int RS_random_seed = 666;
tr1::mt19937 gen(RS_random_seed);
#else
int RS_random_seed = time(0);
tr1::mt19937 gen(RS_random_seed);
#endif // RSDEBUG 
#else
// for some unknown reason, 32-bit compile fails if RS_random_seed is passed as parameter (as above)...
#if RSDEBUG
int RS_random_seed = 666;
tr1::mt19937 gen(666);
#else
// ...there is thus the extremely low possibiity that the value in RS_random_seed and
// the value used to seed the random number stream may differ by 1
int RS_random_seed = time(0);
tr1::mt19937 gen(time(0));
#endif // RSDEBUG 
#endif // RSWIN64 
#endif // !LINUX_CLUSTER

RSrandom::RSrandom(void)
{

#if BATCH && RSDEBUG
DEBUGLOG << "RSrandom::RSrandom(): RS_random_seed=" << RS_random_seed
	<< endl;
#endif // RSDEBUG 

#if RSDEBUG
// RS random initialisation log added by SCFP 25/8/16
//string name = paramsSim->getDir(2) + "RandomLog.txt";
//RSRANDOMLOG.open(name.c_str());
//RSRANDOMLOG << "RSrandom::RSrandom(): creating new random stream" << endl;
#endif
/*
#  if defined(BOOST_HAS_CPP_0X)
cout << "BOOST_HAS_CPP_0X" << endl;
#  else
cout << "***NOT*** BOOST_HAS_CPP_0X" << endl;
#  endif
#  if defined(BOOST_HAS_TR1_RANDOM)
cout << "BOOST_HAS_TR1_RANDOM" << endl;
#     ifdef BOOST_HAS_INCLUDE_NEXT
cout << "BOOST_HAS_INCLUDE_NEXT" << endl;
//#        include_next <random>
#     else
cout << "***NOT*** BOOST_HAS_INCLUDE_NEXT" << endl;
cout << "BOOST_TR1_STD_HEADER" << endl;
//#        include BOOST_TR1_STD_HEADER(random)
#     endif
#  else
cout << "***NOT*** BOOST_HAS_TR1_RANDOM" << endl;
#  endif
*/

normal_x2_valid = 0;
#if !LINUX_CLUSTER
// Set up standard uniform distribution
pRandom01 = new	tr1::uniform_real<> (0.0,1.0);
// Set up standard normal distribution
pNormal = new	tr1::normal_distribution<> (0.0,1.0);
#endif

//cout << endl;
//for (int i = 0; i < 5; i++) {
//	cout << gen() << " ";
//}
//cout << endl << endl;

#if RSDEBUG
//RSRANDOMLOG.close(); RSRANDOMLOG.clear();
#endif
}

RSrandom::~RSrandom(void) {
#if !LINUX_CLUSTER
if (pRandom01 != 0) delete pRandom01;
if (pNormal != 0) delete pNormal;
#endif
}

double RSrandom::Random(void) {
#if LINUX_CLUSTER
return unif_rand();
#else
return pRandom01->operator()(gen);
//double result_Random;
//result_Random = pRandom01->operator()(*gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Random(): result_Random=" << result_Random
//	<< endl;
//#endif
//return result_Random;
#endif
}

int RSrandom::IRandom(int min,int max) {
//#if LINUX_CLUSTER
//return irand(min,max);
//#else
//tr1::uniform_int<> dist(min,max);
//return dist(gen);
//#endif
// output random integer in the interval min <= x <= max (copied from mersenne.cpp)
int r;
r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
if (r > max) r = max;
if (max < min) return 0x80000000;
return r;
}

int RSrandom::Bernoulli(double p) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Bernoulli(): p=" << p
//	<< endl;
#endif
#if LINUX_CLUSTER
return unif_rand() < p;
#else
return Random() < p;
#endif
}

double RSrandom::Normal(double mean,double sd) {
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): mean=" << mean << " sd=" << sd
//	<< endl;
#endif
#if LINUX_CLUSTER
// normal distribution derived from Agner Fog's method
double normal_x1;          // first random coordinate (normal_x2 is member of class)
double w;                  // radius
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x2_valid=" << normal_x2_valid
//	<< endl;
#endif
if (normal_x2_valid) {     // we have a valid result from last call
	normal_x2_valid = 0;
	return normal_x2 * sd + mean;
}
// make two normally distributed variates by Box-Muller transformation
do {
	normal_x1 = 2. * Random() - 1.;
	normal_x2 = 2. * Random() - 1.;
	w = normal_x1*normal_x1 + normal_x2*normal_x2;
} while (w >= 1. || w < 1E-30);
w = sqrt(log(w)*(-2./w));
normal_x1 *= w;  normal_x2 *= w;     // normal_x1 and normal_x2 are independent normally distributed variates
normal_x2_valid = 1;                 // save normal_x2 for next call
#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): normal_x1=" << normal_x1
//	<< endl;
#endif
return normal_x1 * sd + mean;
#else
//double norm = pNormal->operator()(gen);
//#if RSDEBUG
//DEBUGLOG << "RSrandom::Normal(): norm=" << norm
//	<< endl;
//#endif
//return mean + sd * norm;
return mean + sd * pNormal->operator()(gen);
#endif
}

int RSrandom::Poisson(double mean) {
#if LINUX_CLUSTER
return rpois(mean);
#else
tr1::poisson_distribution<> poiss(mean);
return poiss(gen);
#endif
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
