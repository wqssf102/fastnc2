# fastnc
&emsp;&emsp;根据文献报道，从网络图中提取邻阶矩阵（0、1矩阵，0代表指标间没有边，1代表指标间有边），然后依次删除节点并求剩余矩阵的“特征值”。通过将计算结果与删除节点比例
做线性回归并提取slope值，比较不同网络图之间的slope数值大小（简单说，就是看回归线下降的趋势）以反映网络的抗毁性。在实践中，一个网络图有N个节点，那么当
删除M个节点时，M个节点有多个情况的节点组合，目前一些文献是通过degree、BC拓扑特征对节点排序，然后依次删除。在这里，我们用另外的做法，即随机抽取M个节点并
删除、同时删掉剩余矩阵中孤立的节点，正如前面所言，M个节点有多种情况组合，因此我们提供了一个参数n，用来多次迭代删除，如我们在删除123个节点，那么我们从总
节点中随机删除123个节点、删除剩余矩阵孤立节点并求“特征值”，这一步骤重复n次(如1000)，最后我们将n次的结算结果取均值，这样得到的结果更具代表性。  
&emsp;&emsp;在实际数据分析过程中，用随机删除M个节点并迭代n次的做法具有很大的计算量。一开始我们分别用R并行计算、Python的numpy模块计算（同时也用了并行），在给出
指定的线程里计算，结果服务器显示任务超线程了（如给出28线程，结果服务器显示已经调用了500多线程或800多）。然而在不并行的情况下，用R计算300个节点的网络图
的抗毁性需要几小时左右（当然取决于计算机的性能，这里不是说R和Python不好）。最后，我们转向了C++软件，调用了[Eigen库](https://eigen.tuxfamily.org/index.php?title=Main_Page)、
通过[OpenMP](https://www.openmp.org/)并行和通过intel公司的[mkl库](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.esqumy)来加速，300个节点网络的计算在28线程、ram为64G的服务器单节点上，只需要3分钟左右即可算完（迭代1000次）。  

&emsp;&emsp;特别说明，在将C++代码封装为命令行结构时，我们参考了[fastspar](https://github.com/scwatts/fastspar)软件（用C++实现sparcc算法的快速软件）。软件可通过conda安装，
```
conda create -n fastnc -c wqssf102 -c conda-forge -c intel fastnc -y
```
软件的参数如下：
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
  &emsp;&emsp;参数进一步说明：
 ```
-c 为邻接矩阵，必须给出

-o 为输出文件，必须给出

-t 为删除物种（节点的阈值），默认为0.8，即从1%删到80%

-s 为步数，默认为1，假如step为3，那么就是依次删除1个节点，4个节点，7个节点，......,即第N次是删除（N+step）个节点（N>1）
因此，当step大于1时，计算时间将大大缩短，适用于网络图节点多的情况

-g 删除边的比例，即使，先删掉节点，在剩余矩阵中随机删除边

-n 为迭代次数，默认1000，因为我们不知道环境变化会使哪些物种消失，因此这里采用了随机抽样，建议至少迭代1000次。

-j 线程数，默认为4。
 ```
 &emsp;&emsp;邻接矩阵如何获得，方法之一是通过R软件的igraph包导出，如：
 ```
library(igraph)
##gg为igraph对象的网络图数据
adjtab <- as.matrix(get.adjacency(gg))
####导出01矩阵，然后使用服务器计算
write.table(adjtab,"mydir/adj_tab.txt",quote=F,row.names=F,col.names=F,sep="\t")
 ```
  运行如下：
 ```
 conda activate fastnc
time fastnc -c filedir/adj_tab.txt -j 28 -n 1000 -o outdir/ncres.txt
conda deactivate
 ``` 
 &emsp;&emsp;最终结果如下：
 ```
 nc_index	del_numbel
0.579490466952	0.002941176471
0.575350980163	0.005882352941
0.571343702435	0.008823529412
0.568237121284	0.011764705882
0.565799163938	0.014705882353
0.563750376701	0.017647058824
0.562035103083	0.020588235294
0.560567979574	0.023529411765
0.559203409255	0.026470588235
0.557888396978	0.029411764706
0.556521071494	0.032352941176
0.555493312657	0.035294117647
0.554323430896	0.038235294118
 ```
nc_index为删除M个节点后剩余矩阵的“特征值”，del_number为删除节点所占总节点的百分比。  
 
 &emsp;&emsp;最后用R绘图：
 ```
library(ggpmisc)
grp1 <- read.csv("ck.txt",header = T,sep="\t")
grp1$grp <- "CK"
grp2 <- read.csv("NP.txt",header = T,sep="\t")
grp2$grp <- "NP"

####将每组合并
grpnc <- rbind(grp1,grp2)
##指定图例的顺序
grpnc$grp <- factor(grpnc$grp,levels = c("CK","NP"))
 ggplot(grpnc, aes(del_number, nc_index,color=grp)) +##grp为将多个网络图合并时的分组
  # geom_point() +
  geom_smooth(formula = y~x,se = FALSE,method = "lm",show.legend = F)+
  stat_poly_eq(formula = y~x,size=5,family="serif",method = "lm",
               output.type = "numeric",
               parse = T,label.x = 0.3,hjust=0,
               mapping =aes(label = paste(paste(sprintf("italic(R)^2~`=`~%.2g",after_stat(`r.squared`)),"~",
                                                ifelse(after_stat(`p.value`)<=0.001,"'***'",
                                                       ifelse(after_stat(`p.value`)<=0.01,"'**'",
                                                              ifelse(after_stat(`p.value`)<=0.05,"'*'",NA)))),
                                          sprintf( "italic(Slope)~`=`~%.3g",after_stat(b_1)),sep = "*\", \"*")))+
  labs(x="Proportion of removed nodes",
       y="Natural connectivity",
       color="Group")+
  scale_color_manual(values = c("CK"="#7570B3","NP"="#1B9E77"))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 20,colour = "black"))+
  theme(axis.text = element_text(colour = "black",size = 20,angle = 0,vjust = 0.5))+
  theme(legend.title =  element_text(size = 20,colour = "black"), 
        legend.text = element_text(size = 20,colour = "black"),
        strip.text =  element_text(size = 20,colour = "black"),
        text = element_text(family = "serif"))


 ```
 ## 关于计算速度和时间
 &emsp;&emsp; 1、我们建议将网络图的节点控制在500个以下，网络节点太多导致计算量太大，服务器算不过来。
 &emsp;&emsp; 2、若用户的网络包含的节点比较多（如几千个节点），假如用户想要迭代1000次，当服务器的运行内存不大时，我们建议分为4个任务提交到服务器中，每个任务执行250次的迭代，最后取4个任务计算结果的均值即可，这样会大大缩减用户等待的时间。特别地，可以通过设置step参数（如step设置为3）来扩大删除节点的步长，这样耗时将显著减少。总之，step是一个很有用的参数。我们的软件还在开发中，若你对软件的功能有需求或发现软件问题，请联系作者：565715597@qq.com  
 &emsp;&emsp;fastnc软件目前只在Linux系统下测试，其他系统没测试过（Win系统肯定是不行）。src为软件的源码，源码需要编译才能使用，用户若想自己编译，那么需要配置依赖的库和软件,编译方式
为：
```
MKLROOT=/opt/intel/oneapi/mkl/2022.1.0
##
g++ -std=c++11 -O3 -fopenmp -march=native -mavx -mfma -o fastnc fastnc.cpp fastnc_opts.cpp common.cpp -DMKL_ILP64 -m64\
 ${MKLROOT}/lib/intel64/libmkl_scalapack_ilp64.a \
 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_cdft_core.a \
 ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a \
 ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a \
 ${MKLROOT}/lib/intel64/libmkl_core.a\
 ${MKLROOT}/lib/intel64/libmkl_blacs_openmpi_ilp64.a \
 -Wl,--end-group -lgomp -lpthread -lm -ldl
``` 

 
 参考文献：  
 
 [1] Wu MH, Chen SY, Chen JW, Xue K, Chen SL, Wang XM, Chen T, Kang SC, Rui JP, Thies JE, Bardgett RD, Wang YF. Reduced microbial stability in the active layer is associated with carbon loss under alpine permafrost degradation. Proc Natl Acad Sci U S A. 2021 Jun 22;118(25):e2025321118. doi: 10.1073/pnas.2025321118.  
 
[2] G.-s. Peng, J. Wu, Optimal network topology for structural robustness based on natural connectivity. Physica A 443, 212–220 (2016). doi: 10.1016/j.physa.2015.09.023
 
 
 
 
  
  
