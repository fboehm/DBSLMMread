/*
Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)
Copyright (C) 2019  Sheng Yang and Xiang Zhou

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <math.h>

#include "dtpr.hpp"
#include "dbslmmfit.hpp"
#include "dbslmm.hpp"
#include "tobool.h"


#include <armadillo>

#include "../include/dtpr.hpp"
#include "../include/dbslmmfit.hpp"
#include "../include/dbslmm.hpp"

using namespace std;
using namespace arma;

DBSLMM::DBSLMM(void) :
	version("0.3"), date("05/01/2021"), year("2021")
{}

void DBSLMM::printHeader(void)
{
	cout << endl;
	cout << "*************************************************************"<< endl;
	cout << "  Deterministic Bayesian Sparse Linear Mixed Model (DBSLMM)  " << endl;
	cout << "  Version " << version << ", " << date << "                  " << endl;
	cout << "  Visit http://www.xzlab.org/software.html For Update        " << endl;
	cout << "  (C) " << year << " Sheng Yang, Xiang Zhou                  " << endl;
	cout << "  GNU General Public License                                 " << endl;
	cout << "  For Help, Type ./dbslmm -h                                 " << endl;
	cout << "*************************************************************" << endl;
	cout << endl;

	return;
}

void DBSLMM::printHelp(void) {
	cout << " FILE I/O RELATED OPTIONS" << endl;
	cout << " -s        [filename]  " << " specify input the summary data for the small effect SNPs." << endl;
	cout << " -l        [filename]  " << " specify input the summary data for the large effect SNPs." << endl;
	cout << " -r        [filename]  " << " specify input the bfile of reference data." << endl;
	cout << " -n        [num]       " << " specify input the sample size of the summary data." << endl;
	cout << " -mafMax   [num]       " << " specify input the maximium of the difference between reference panel and summary data." << endl;
	cout << " -nsnp     [num]  " << " specify input the number of snp." << endl;
	cout << " -b        [num]       " << " specify input the block information." << endl;
	cout << " -h        [num]       " << " specify input the heritability." << endl;
	cout << " -t        [filename]  " << " specify input thread." << endl;
	cout << " -eff      [filename]  " << " specify output the estimate effect SNPs." << endl;
	return;
}

void DBSLMM::Assign(int argc, char ** argv, PARAM &cPar) {
	
	string str;
	for (int i = 0; i < argc; i++) {

		if (strcmp(argv[i], "--smallEff") == 0 || strcmp(argv[i], "-s") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.s = str;
		}
		else if (strcmp(argv[i], "--largeEff") == 0 || strcmp(argv[i], "-l") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.l = str;
		}
		else if (strcmp(argv[i], "--reference") == 0 || strcmp(argv[i], "-r") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.r = str;
		}
		else if (strcmp(argv[i], "--N") == 0 || strcmp(argv[i], "-n") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.n = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--mafMax") == 0 || strcmp(argv[i], "-mafMax") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.mafMax = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--numSNP") == 0 || strcmp(argv[i], "-nsnp") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.nsnp = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--block") == 0 || strcmp(argv[i], "-b") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.b = str;
		}
		else if (strcmp(argv[i], "--Heritability") == 0 || strcmp(argv[i], "-h") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.h = atof(str.c_str());
		}
		else if (strcmp(argv[i], "--Thread") == 0 || strcmp(argv[i], "-t") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.t = atoi(str.c_str());
		}
		else if (strcmp(argv[i], "--EFF") == 0 || strcmp(argv[i], "-eff") == 0) {

			if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
			++i;
			str.clear();
			str.assign(argv[i]);
			cPar.eff = str;
		}
		else if (strcmp(argv[i], "--training") == 0 || strcmp(argv[i], "-training") == 0) {
		  
		  if (argv[i + 1] == NULL || argv[i + 1][0] == '-') { continue; }
		  ++i;
		  str.clear();
		  str.assign(argv[i]);
		  cPar.training = to_bool(str);
		}
		
	}
	return;
}

