#include "fastnc.h"

int main(int argc, char **argv) {
	FastNCOptions fastnc_options = get_commandline_arguments(argc, argv);
	// // Initialise a FastNC object
	FastNC fastnc(fastnc_options.adj_table_filename, fastnc_options.outfile_filename, fastnc_options.threshold_name, fastnc_options.number_name, fastnc_options.jobs_name);
	
	//将natural_connectivity计算1000次，然后矩阵相加，最后除以1000即可
	MatrixXd data_M = readmyfile(fastnc_options.adj_table_filename);
	// MatrixXd data_R = natural_connectivity(data_M, fastnc_options.threshold_name);//先计算依次，得到一个初始矩阵
	// MatrixXd data_T = Matrix<double, Dynamic, 2>();;
	// data_T.resize(ceil(data_M.rows()*fastnc_options.threshold_name), 2);
	// //
	MatrixXd data_R = Matrix<double, Dynamic, 2>();
	data_R.setZero(ceil(data_M.rows()*fastnc_options.threshold_name),2);
		omp_set_num_threads(fastnc_options.jobs_name);
		#pragma omp declare reduction( + : Eigen::MatrixXd : omp_out += omp_in ) initializer( omp_priv = omp_orig )
		#pragma omp parallel for reduction(+:data_R), schedule(static, 1)
	for (int i =1;i <= fastnc_options.number_name;i++)
	{
		printf("\tRunning %d\n", i);
		// MatrixXd data_T = natural_connectivity(data_M, fastnc_options.threshold_name);
		// data_R = data_R + data_T;//将结果相加
		data_R = data_R + natural_connectivity(data_M, fastnc_options.threshold_name);
	}
	// MPI::Finalize();
	 data_R = data_R / fastnc_options.number_name * 1.0;

	
	// for(int i =0; i < data_R.rows(); i++){
		// data_R(i,0) = data_R(i,0) / data_R(i,2);
		// data_R(i,1) = data_R(i,1) / data_R(i,2);
	// }
	ofstream datafile;
	datafile.open(fastnc_options.outfile_filename, ios::out);
	// datafile <<"nc_index"<< "\t" << "del_numbel" << "\t" << "NOTnan_num";
	datafile <<"nc_index"<< "\t" << "del_numbel";
	datafile << std::endl;
	    for (unsigned int i = 0; i < data_R.rows(); ++i) {
        for (unsigned int j = 0; j < 2; ++j) {
			if(j == 0 ){
				datafile << std::fixed << std::setprecision(12) << data_R(i, j)<<"\t";
			}else{
				datafile << std::fixed << std::setprecision(12) << data_R(i, j);
			}
        }
        datafile << std::endl;
    }	
    return 0;
}


FastNC::FastNC(const std::string &inputfile, const std::string &outfile, const float &_threshold, const int &_number, const int &_jobs) {
    adj_table_filename = inputfile;
    outpfile_filename = outfile;
	threshold_name =  _threshold;
	number_name = _number;
	jobs_name = _jobs;
}

////定义一个获取文件行数的函数，用来设置Matrix的行和列
int GetFileCount(const std::string& szFile)
{
	fstream fin(szFile, ios::in);
	if (!fin)
	{
		cerr << "can not open file" << endl;
		return -1;
	}
	char c;
	int lineCnt = 0;
	while (fin.get(c))
	{
 		if (c == '\n')
			lineCnt++;
	}
	//cout << lineCnt << endl;
	fin.close();

	return lineCnt;
}
//
////定义一个拆分文件的函数
void _split(const string& s, char delim, vector<string>& elems) {
	stringstream ss(s);
	string item;

	while (getline(ss, item, delim)) {
		elems.push_back(item);
	}
}
vector<string> split(const string& s, char delim) {
	vector<string> elems;
	_split(s, delim, elems);
	return elems;
}

//读取数据
MatrixXd readmyfile(const std::string& ldfile)
{
	MatrixXd myfile = Matrix<double, Dynamic, Dynamic, ColMajor>();
	int rowCount = GetFileCount(ldfile);
	myfile.resize(rowCount, rowCount);
	ifstream infile(ldfile);
	if (infile.is_open())
	{
		std::string line;
		int row = 0;
		while (getline(infile, line) && row < rowCount)
		{
			//cout << line << endl;
			std::vector<std::string> arrContexts = split(line, '	');
			int column = 0;
			for (auto value : arrContexts)
			{
				if (column < rowCount)
				{
					myfile(row, column) = atof(value.c_str());
				}
				++column;
			}
			++row;
		}
	}

	return myfile;
}

////定义一个计算特征值和进一步计算最终目标的函数
float get_nc(MatrixXd& M, int const tt_num)
{
	VectorXcd eivals = M.eigenvalues();
	MatrixXd res = eivals.real();
        res = res.array().exp();//获得实数部分并且做exp计算
	double rest = log(res.sum() / res.rows());
	rest = rest / ((tt_num*1.0) - log((tt_num*1.0)));
	return rest;
}



////定义一个随机抽样的函数
set<int> mysample(int const total_num, int const del_num)
{
	set<int> del_num_res;
	// srand(time(0));bug
	std::mt19937 rng(std::random_device{}());
	srand(rng());
	for (int i = 1; i <= del_num; i++)
	{
		del_num_res.insert(rand() % total_num + 1);

	}

	//判断生成的数据的个数是否满足要求
	while (del_num_res.size() < del_num)
	{
		int need_num = del_num - del_num_res.size();
		for (int i = 1; i <= need_num; i++)
		{
			std::mt19937 rng(std::random_device{}());
			srand(rng());
			del_num_res.insert(rand() % total_num + 1);
			if (del_num_res.size() == del_num)
			{
				break;
			}
		}
	}
	return del_num_res;
}



////定义一个删除Martix行和列的函数
void remove_row_col(MatrixXd& M, unsigned int rowToRemove) {
	unsigned int numRows = M.rows() - 1;
	unsigned int numCols = M.cols() - 1;

	if (rowToRemove < numRows) {
		M.block(rowToRemove, rowToRemove, numRows - rowToRemove, numCols - rowToRemove) =
			M.block(rowToRemove + 1, rowToRemove + 1, numRows - rowToRemove, numCols - rowToRemove);
	}

	M.conservativeResize(numRows, numCols);
}

//获取行列为0
int getFilterRanks(MatrixXd& matrix, std::set<int>& setRanks)
{
	if (matrix.rows() <= 0 || matrix.cols() <= 0)
		return 0;

	for (int rank = 0; rank < matrix.rows(); ++rank)
	{
		if (matrix.row(rank).sum() == 0)
		{
			setRanks.insert(rank);
		}
		else if (matrix.col(rank).sum() == 0)
		{
			setRanks.insert(rank);
		}
	}
	return 0;
}


bool deleteRanks(MatrixXd& matrix, const std::set<int>& setRanks)
{
	if (matrix.rows() <= 0 || matrix.cols() <= 0)
		return false;
	if (setRanks.size() <= 0)
		return false;
	auto it = setRanks.rbegin();
	while (it != setRanks.rend())
	{
		//std::cout << *it << std::endl;
		remove_row_col(matrix, *it);		// 删除行列
		++it;
	}
	return true;
}

//
MatrixXd natural_connectivity(const MatrixXd& re_matrix,  const float del_per)
{
	int total_nodes = re_matrix.rows();//获取原始矩阵的行数
	MatrixXd resdt = Matrix<double, Dynamic, 2>();//定义一个 total_nodes行、2列的矩阵，用来存放最后计算结果
	resdt.setZero(ceil(total_nodes * del_per), 2);//对定义的矩阵初始化
	//2 获取随机数
	for (int rmnd = 0; rmnd <= (ceil(total_nodes * del_per) - 1); rmnd++)
	{
		//printf("Running:%d\n",rmnd+1);
		MatrixXd matrix = re_matrix;
		set<int> rows = mysample(total_nodes, rmnd + 1);
		//cout << rows.size() << endl;
		//printSet(rows);
		// 
		//3 删除行、列
		bool bRet = deleteRanks(matrix, rows);
		if (bRet)
		{
			//std::cout << "delete successed!" << std::endl;
		}


		//4 过滤剩余行列和为0的
		rows.clear();
		int counts = getFilterRanks(matrix, rows);
		if (counts > 0)
		{
			bRet = deleteRanks(matrix, rows);
			if (bRet)
			{
				//std::cout << "delete successed!" << std::endl;
			}
		}
		//
		// float ncres = get_nc(matrix, total_nodes);
		// float isnonb = isnan(ncres);//is number?
		// if(isnonb==1){
			// resdt(rmnd, 0) = 0.0;
			// resdt(rmnd, 2) = 0.0;
		// }else{
			// resdt(rmnd, 0) = ncres;
			// resdt(rmnd, 2) = 1.0;
		// }
		resdt(rmnd, 0) = get_nc(matrix, total_nodes);
		resdt(rmnd, 1) = (rmnd + 1) * 1.0 / total_nodes * 1.0;
		//resdt(rmnd, 2) = rmnd +1;
	}
	return resdt;
}
