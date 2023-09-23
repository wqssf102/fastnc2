#include "common.h"

std::vector<int> mysample(int const total_num, int const del_num)
{
    std::vector<int> del_num_res;
    std::random_device rd;
    std::mt19937 rng(rd()+omp_get_thread_num());
    std::uniform_int_distribution<int> uni(0, total_num - 1);
    while (del_num_res.size() < del_num)
    {
        int num = uni(rng);
        if (std::find(del_num_res.begin(), del_num_res.end(), num) == del_num_res.end()) {
            del_num_res.push_back(num);
        }
    }
    return del_num_res;
}

std::vector<std::vector<int>> remove_random(std::vector<std::vector<int>>& matrix, std::vector<int>& remove_nodes)
{ 
		std::sort(remove_nodes.begin(), remove_nodes.end(), std::greater<int>());

		for (int i = 0; i < remove_nodes.size(); ++i)
		{
			int index = remove_nodes[i];
			matrix.erase(matrix.begin() + index); // 删除第index行
		}
		for (auto& row : matrix)
		{
			for (int i = 0; i < remove_nodes.size(); ++i)
			{
				int index = remove_nodes[i];
				row.erase(row.begin() + index); // 删除每一行的第index列
			}
		}

	return matrix;
}
// 计算矩阵的自然连通度
double get_nc(std::vector<std::vector<int>>& matrix)
{
    const int n = matrix.size(); // 矩阵维度
    std::vector<double> A(n * n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i * n + j] = matrix[i][j];
        }
    }

    std::vector<double> eigenvalues(n); // 特征值
    MKL_INT info;
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'N', 'U', n, A.data(), n, eigenvalues.data());
	double sum = 0;
	if (info == 0) {
    std::vector<double> expValues(n); // 存储特征值的exp计算结果
    vdExp(n, eigenvalues.data(), expValues.data()); // 对特征值进行exp计算
    sum = cblas_dasum(n, expValues.data(), 1);
    sum = log(sum / n);
    sum = sum / (n - log(n));
	}else{
	std::cout << "Failed to compute eigenvalues." << std::endl;
	}
    return sum;
}
// 组合函数
matrix<double, row_major> natural_connectivity(std::vector<std::vector<int>>& matrix_data, const float del_per)
{
    int total_nodes = matrix_data.size(); // 获取原始矩阵的行数
    int num_del = ceil(total_nodes * del_per); // 计算需要删除的节点数
    matrix<double, row_major> resdt(num_del, 2);
    for (int rmnd = 0; rmnd < num_del; rmnd++)
    {
        std::vector<std::vector<int>> matrix = matrix_data;
        if ((total_nodes - rmnd) < 3)
        {
            printf("The number of nodes is too small\n");
            exit(1);
        }
        // 随机抽节点
        std::vector<int> rows = mysample(total_nodes, rmnd + 1);
        // 删除节点
        std::vector<std::vector<int>> matrix_ct = remove_random(matrix, rows);
        // 计算自然连通度
        resdt(rmnd, 0) = get_nc(matrix_ct);
        resdt(rmnd, 1) = (rmnd + 1) * 1.0 / total_nodes * 1.0;
    }
    return resdt;
}
// 当使用步长
matrix<double, row_major> natural_connectivity_step(std::vector<std::vector<int>>& matrix_data, const float del_per,const int step)
{
	// matrix_data.size为原始矩阵的行数，但后面的代码是从0开始，因此需要-1
	int total_nodes = matrix_data.size();//获取原始矩阵的行数
	int del_total_nodes = ceil(total_nodes * del_per);
	int mtrows = del_total_nodes / step ;
	if(mtrows < 2)
	{
			printf("The number of step is too big\n");
			exit(1);
	}
	
	// if((del_total_nodes-(mtrows * step)) > 0)
	// {
		// mtrows = mtrows + 1;
	// }
	 matrix<double, row_major> resdt(mtrows, 2);

for (int rmnd = 0; rmnd < del_total_nodes; rmnd=rmnd+step)
	{
		std::vector<std::vector<int>> matrix = matrix_data;
		//判断最后一个循环的结果是否大于删除总的节点数，若大于，则直接等于总节点数
		if (rmnd + step > del_total_nodes)
			{
				rmnd = del_total_nodes ;
			}
		if ((total_nodes - rmnd) < 3)
		{
			printf("The number of nodes is too small\n");
			exit(1);
		}
		 std::vector<int> rows = mysample(total_nodes, rmnd + 1);
		 std::vector<std::vector<int>> matrix_ct = remove_random(matrix, rows);		
		if (rmnd + step > del_total_nodes)
			{
			resdt((resdt.size1())-1, 0) = get_nc(matrix_ct);
		    resdt((resdt.size1())-1, 1) = (rmnd + 1) * 1.0 / total_nodes * 1.0;
			}else{
			resdt(rmnd/step, 0) = get_nc(matrix_ct);
		    resdt(rmnd/step, 1) = (rmnd + 1) * 1.0 / total_nodes * 1.0;
			}

	}
	return resdt;
}

std::vector<std::vector<int>> read_file(const std::string& filename)
{
	std::vector<std::vector<int>> matrix_data;
    std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line)) {
            std::vector<int> row;
            std::istringstream iss(line);
            int value;
            while (iss >> value) {
                row.push_back(value);
            }
            matrix_data.push_back(row);
        }
        file.close();
	return matrix_data;
}

void write_tsv(const std::string& outfile, matrix<double> data_R) {
	ofstream datafile;
	datafile.open(outfile, ios::out);
	datafile <<"nc_index"<< "\t" << "del_number" << "\n";
		for (size_t i = 0; i < data_R.size1(); ++i) {
			for (size_t j = 0; j < 2; ++j) {
				datafile << std::setprecision(6) <<std::fixed << data_R(i, j);
				if (j < 1)
					datafile << "\t";
			}
			datafile << "\n";
		}
}