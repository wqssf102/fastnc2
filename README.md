# FastNC2：Fast calculation of natural connectivity
## Installing
FastNC2 can be installed using conda or from source.
## Conda
To install through conda, use:
```
conda create -n fastnc2 -c wqssf102 fastnc2
```
## Or download the compiled version
```
1. Download: https://github.com/wqssf102/fastnc2/releases/download/fastnc2/fastnc2.tar.gz
2. tar xzvf fastnc2.tar.gz
3. chmod a+x fastnc2
```
**help**：
```
fastnc2 --h
Program: FastNC2 (use c++ to calculate the natural connectivity)
Contact: Qiusheng WU (565715597@qq.com)

  -i [ --input ] arg          input.txt,The result from get.adjacency function of igraph package
  -o [ --output ] arg         output.txt
  -s [ --step ] arg (=1)      number of threads (default 1)
  -t [ --thread ] arg (=0.80) the threshold for deletion of node (default: 0.8)
  -n [ --number ] arg (=1000)  number of iterations (default: 1000)
  -j [ --jobs ] arg (=4)       number of jobs (default: 4)
  ```

## USAGE and Tutorial

* Chinese version

  <u>https://wqssf102.github.io/fastnc2/</u>

* English version

  In the plan
  
  
