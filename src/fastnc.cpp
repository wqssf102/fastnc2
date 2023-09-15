#include "common.h"
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

//main函数
int main(int argc, char **argv) {
    // 声明命令行参数选项
	double default_thread = 0.8;
    double thread;
	int default_step = 1;
    int step;
	int default_number = 1000;
    int number;
	int default_jobs = 4;
    int jobs;
    po::options_description desc("  Program: FastNC (use c++ to calculate the natural connectivity).\n" 
		"  Contact: Qiusheng WU (565715597@qq.com)");
    desc.add_options()
        ("help,h", "help message")
        ("input,i", po::value<std::string>(), "input.txt,The result from get.adjacency function of igraph package")
		("output,o", po::value<std::string>(), "output.txt")
		("step,s", po::value<int>(&step)->default_value(default_step), "number of threads (default 1)")
		("thread,t", po::value<double>(&thread)->default_value(default_thread,"0.80"), "the threshold for deletion of node (default: 0.8)")
		("number,n", po::value<int>(&number)->default_value(default_number), " number of iterations (default: 1000)")
		("jobs,j", po::value<int>(&jobs)->default_value(default_jobs), " number of jobs (default: 4)");

				
    // 解析命令行参数选项
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 0;
    }
    if (vm.empty()) {
        std::cout << desc << "\n";
        return 0;
    }
	//
	int coreNum = omp_get_num_procs();//获得处理器个数
	if(coreNum > 1 &&  jobs > coreNum){
		fprintf(stderr, "\n%s: error: only %d threads are available\n",  argv[0], coreNum); 
		exit(1);
	}
	// 读取数据
	std::vector<std::vector<int>> matrix_data = read_file(vm["input"].as<string>());
if(vm["step"].as<int>()>1){
	int mtrows_tt_tmp = ceil(matrix_data.size()*vm["thread"].as<double>());
	int mtrows_tt = mtrows_tt_tmp / vm["step"].as<int>();
	// if((mtrows_tt_tmp-(mtrows_tt * vm["thread"].as<double>())) > 0)
	// {
		// mtrows_tt = mtrows_tt + 1;
	// }
	matrix<double> data_R(mtrows_tt, 2);
	omp_set_num_threads(vm["jobs"].as<int>());
	#pragma omp declare reduction( + : boost::numeric::ublas::matrix<double> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
    #pragma omp parallel for reduction(+: data_R)
	for (int i =1;i <= vm["number"].as<int>();i++)
	{
		printf("\tRunning %d\n", i);
		data_R = data_R + natural_connectivity_step(matrix_data, vm["thread"].as<double>(),vm["step"].as<int>());
	}
	
	data_R = data_R / vm["number"].as<int>() * 1.0;//取跑N次的均值
		// 输出结果
	write_tsv(vm["output"].as<std::string>(),data_R);

}else{
	matrix<double> data_R(ceil(matrix_data.size() * vm["thread"].as<double>()), 2);
	omp_set_num_threads(vm["jobs"].as<int>());
	#pragma omp declare reduction( + : boost::numeric::ublas::matrix<double> : omp_out += omp_in ) initializer( omp_priv = omp_orig )
    #pragma omp parallel for reduction(+: data_R)
	for (int i =1;i <= vm["number"].as<int>();i++)
	{
		printf("\tRunning %d\n", i);
		data_R = data_R + natural_connectivity(matrix_data, vm["thread"].as<double>());
	}
	
	data_R = data_R / vm["number"].as<int>() * 1.0;//取跑N次的均值
		// 输出结果
	write_tsv(vm["output"].as<std::string>(),data_R);
}



	return 0;
}	
