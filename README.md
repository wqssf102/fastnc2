# fastnc：Fast calculation of natural connectivity
## INSTALLATION
**conda**
```
conda create -n fastnc -c wqssf102 -c conda-forge -c intel fastnc -y
```
**help**：
```
fastnc --h
Program: FastNC (use c++ to calculate the natural connectivity)
Contact: Qiusheng WU (565715597@qq.com)

  fastnc [options] --adj_table <file> --outfile <file>
  -c <file>, --adj_table <file>
                The result from get.adjacency function of igraph package
  -o <file>, -outfile <file>
                Result output table
Options:
  -t <float>, -threshold <float>
                The threshold for deletion of node (default: 0.8)
  -s <int>, -step <int>
                The step for deletion of node (default: 1)
  -g <float>, -edge <float>
                The threshold for deletion of edge (default: 0)
  -n <int>, -number <int>
                Number of iterations (default: 1000)
  -j <int>, -job <int>
                Number of jobs (default: 4)
Other:
  -h        --help
                Display this help and exit
  ```

## USAGE and Tutorial

* Chinese version

  <u>https://wqssf102.github.io/fastnc/</u>

* English version

  In the plan
  
  
