#ifndef __FASTNC_H__
#define __FASTNC_H__

#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif

#ifndef EIGEN_VECTORIZE_SSE4_2
#define EIGEN_VECTORIZE_SSE4_2
#endif

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <stdio.h> 
#include <cstdio>
#include <random>
#include <vector>
#include <algorithm>
#include <math.h>
#include <set>
#include <getopt.h>
#include <cstdlib>
#include <ctime>
#include "omp.h"
#include "common.h"

#include "fastnc_opts.h"



using namespace Eigen;
using namespace std;

//定义一个名为FastNC的结构体
//结构体FastNC
class FastNC {
public:
	FastNC(const std::string &inputfile, const std::string &outfile, const float &_threshold, const int &_number, const int &_jobs);
	MatrixXd natural_connectivity();
private:
	std::string adj_table_filename;
    std::string outpfile_filename;
	float threshold_name;
	int number_name;
	int jobs_name;
	};


//

int GetFileCount(const std::string& szFile);
void _split(const string& s, char delim, vector<string>& elems);
vector<string> split(const string& s, char delim);
MatrixXd readmyfile(const std::string& inputfile);
set<int> mysample(int const total_num, int const del_num);
void remove_row_col(MatrixXd& M, unsigned int rowToRemove);
int getFilterRanks(MatrixXd& matrix, std::set<int>& setRanks);
bool deleteRanks(MatrixXd& matrix, const std::set<int>& setRanks);
float get_nc(MatrixXd& M, int const tt_num);
MatrixXd natural_connectivity(const MatrixXd& re_matrix,  const float del_per);

#endif