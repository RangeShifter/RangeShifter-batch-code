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
 
 
/*------------------------------------------------------------------------------

RangeShifter v2.0 RSrandom

Implements the RSrandom class

Author: Steve Palmer, University of Aberdeen

Last updated: 9 November 2020 by Steve Palmer

------------------------------------------------------------------------------*/

#ifndef RSrandomH
#define RSrandomH

#include <stdlib.h>
#include <fstream>
//#include <iostream>

#include "Version.h"

using namespace std;

#if RSDEBUG
extern ofstream DEBUGLOG;
#endif

	#if LINUX_CLUSTER
	//#include <random>
	//#include <tr1/random>
	#include "maths.h"
	#else
	#if RSWIN64
	#include <dinkumware64/random>
	#else
	#include <dinkumware/random>
	#endif
	#endif

	class RSrandom {

	public:
		RSrandom(void);
		~RSrandom(void);
		double Random(void);
		int IRandom(int,int);
		int Bernoulli(double);
		double Normal(double,double);
		int Poisson(double);

	private:
		double normal_x2; int normal_x2_valid; // variables used by Normal distribution
	#if !LINUX_CLUSTER
		tr1::uniform_real<> *pRandom01;
		tr1::normal_distribution<> *pNormal;
	#endif
	};

//---------------------------------------------------------------------------
#endif // RSrandomH
