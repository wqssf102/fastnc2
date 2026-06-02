#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <omp.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <mkl.h>
#include <mkl_lapacke.h>
#include <mkl_vml.h>
#include <mkl_cblas.h>

using namespace boost::numeric::ublas;
using namespace std;
std::vector<int> mysample(int const total_num, int const del_num);
std::vector<std::vector<int>> remove_random(std::vector<std::vector<int>>& matrix, std::vector<int>& remove_nodes);
double get_nc(std::vector<std::vector<int>>& matrix);
matrix<double, row_major> natural_connectivity(std::vector<std::vector<int>>& matrix_data, const float del_per);
matrix<double, row_major> natural_connectivity_step(std::vector<std::vector<int>>& matrix_data, const float del_per,const int step);
std::vector<std::vector<int>> read_file(const std::string& filename);
void write_tsv(const std::string& outfile, matrix<double> data_R);
