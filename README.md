# fastnc
&emsp;&emsp;从网络图的邻阶矩阵（0、1矩阵，0代表指标间没有边，1代表指标间有边）中删除节点并求剩余矩阵“特征值”被用来表征网络抗毁性，该方法被用在比较不同组的微生物共现网络图的稳定性中。用法为：删除网络图中的1、2、3、4、5、6、7、......N个节点并求取每次删后矩阵的特征值并做转化，将最终结果称为自然连通度。在实际中，一个网络图包含多个节点，那么如何使得某次删除的节点具有代表性呢？一种方法是根据节点的Degree大小来排序，Degree从大到小依次删除节点并计算剩余节点的邻接矩阵的自然连通度。然而，当网络图中的大部分节点的Degree比较接近甚至相等时，按照Degree排序删除节点的做法将不再有效。实施上，我们可以近似认为用删除节点的方式是模拟环境变化使得群落中某些物种消失（或某些物种之间的关系消失，即删除边）之后评估群落的稳定性如何。那么，我们事先不知道哪些物种会消失，例如，当群落中的N个物种消失时，我们很难知道这些N个物种的组合里具体是哪些物种，因此我们使用随机抽样的方法对节点进行随机抓取并删除的方法来模拟物种的随机消失，即计算随机抓取N个节点并删除之后的自然连通度。为了提高结果的可重复性，我们设定删除N个节点时可以迭代1000次（甚至更多次）并分别计算每次删除之后的自然连通度，最后取均值。当网络图节点数目大时，邻接矩阵特征值的计算量非常大，我们在实际计算过程中发现，当删除N个节点时，其自然连通度与（N-1）和（N+1）的结果很接近，因此，在尽量不改变原始结果的情况下，我们设定了步长（step）。当步长大于1，假设为3，那么删除的节点数将为1、4、7、10......（N+step），在这种设置下，计算所需时间将比原始数据（步长为1）至少减少3倍。  
&emsp;&emsp;基于上述内容，我们基于C++语言和调用了[Eigen库](https://eigen.tuxfamily.org/index.php?title=Main_Page)、
通过[OpenMP](https://www.openmp.org/)并行和[mkl库](https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html#gs.esqumy)来加速，300个节点网络的计算在28线程、ram为64G的服务器单节点上，只需要3分钟左右即可算完（迭代1000次）。  

&emsp;&emsp;特别说明，在将C++代码封装为命令行结构时，我们参考了[fastspar](https://github.com/scwatts/fastspar)软件（用C++实现sparcc算法的快速软件）。fastnc可通过conda安装，
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
setwd("D:\\fastnc")
##读取物种数据
sp <- read.csv("sp.csv",header = T,row.names = 1)
##sp的数据结构
head(sp,c(2,2))
    g__Candidatus.Udaeobacter g__uncultured__f__uncultured__o__Acidobacteriales
site1               0.001288199                                       0.001502898
site2               0.002008320                                       0.000788983
##读取分组文件
grp <- read.csv("grp.csv",header = T,row.names = 1)
##样本分组文件
head(grp)
      sampleid group
site1    site1    AA
site2    site2    AA
site3    site3    AA
site4    site4    AA
site5    site5    AA
site6    site6    AA
##
sp <- sp[grp$sampleid,]
#######按照分组计算物种之间的关系
library(dplyr)
library(igraph)
##calcor.R为已经整理好的代码，方便计算相关性
source("calcor.R",encoding = "utf-8")
corres <- calcor(codt = sp,group = grp,r_th = 0.6,p_th = 0.5,type = "cor")
##提取网络图的邻接矩阵
for (i in unique(grp$group)) {
  adjtab <- as.matrix(get.adjacency(corres[[i]]$gg))
  ####导出01矩阵，然后使用服务器计算
  write.table(adjtab,paste("NC/",i,".txt",sep = ""),quote=F,row.names=F,col.names=F,sep="\t")
}
 ```
  运行如下：
 ```
##激活环境
conda activate fastnc 
##开始计算，当网络总结点数除以step（即-s参数）小于100时候。不建议使用step，默认即可。
##因为当step大于1时，软件会按照步长删除节点，步长太大，最后剩余的节点过少，就会影响最终绘图时拟合的规律
time fastnc -c AA.txt -j 28 -n 1000 -s 3 -o AA_res.txt 
#real   2m25.998s
#user   64m28.868s
#sys    0m7.128s
##
time fastnc -c BB.txt -j 28 -n 1000 -s 3 -o BB_res.txt 
#real   2m30.020s
#user   65m51.200s
#sys    0m8.165s
##关闭conda环境
conda deactivate
##在实际跑数据中，step是否使用取决于服务器的配置，其次是节点的数量。
 ``` 
 &emsp;&emsp;最终结果如下：
 ```
nc_index    del_number
0.376387425101  0.002493765586
0.368335959568  0.009975062344
0.361530920709  0.017456359102
0.356399175783  0.024937655860
0.352021378405  0.032418952618
 ```
nc_index为删除M个节点后剩余矩阵的“特征值”，del_number为删除节点所占总节点的百分比。  
 
 &emsp;&emsp;最后用R绘图：
 ```
library(ggpmisc)
grp1 <- read.csv("NC/AA_res.txt",header = T,sep="\t")
grp1$grp <- "AA"
grp2 <- read.csv("NC/BB_res.txt",header = T,sep="\t")
grp2$grp <- "BB"
####将每组合并
grpnc <- rbind(grp1,grp2)
##指定图例的顺序
grpnc$grp <- factor(grpnc$grp,levels = c("AA","BB"))
ggplot(grpnc, aes(del_number, nc_index,color=grp)) +##grp为将多个网络图合并时的分组
  geom_point() +
  geom_smooth(formula = y~x,se = FALSE,method = "lm",show.legend = F)+
  stat_poly_eq(formula = y~x,size=5,family="serif",method = "lm",
               output.type = "numeric",
               parse = T,label.x = 0.3,hjust=0,
               mapping =aes(label = paste("italic(R)^2~`=`","~",sprintf("\"%#.*f\"",3,after_stat(`r.squared`)),"*\"\"*",
                ifelse(after_stat(`p.value`)<=0.001,"'***'",
                ifelse(after_stat(`p.value`)<=0.01,"'**'",
                 ifelse(after_stat(`p.value`)<=0.05,"'*'",NA))),"*\", \"*",
                  "italic(Slope)~`=`","~",sprintf("\"%#.*f\"",3,after_stat(b_1)),sep = "")))+
  labs(x="Proportion of removed nodes",
       y="Natural connectivity",
       color="Group")+
  scale_color_manual(values = c("AA"="#7570B3","BB"="#1B9E77"))+
  theme_bw()+theme(panel.grid = element_blank())+
  theme(axis.title = element_text(size = 20,colour = "black"))+
  theme(axis.text = element_text(colour = "black",size = 20,angle = 0,vjust = 0.5))+
  theme(legend.title =  element_text(size = 20,colour = "black"), 
        legend.text = element_text(size = 20,colour = "black"),
        strip.text =  element_text(size = 20,colour = "black"),
        text = element_text(family = "serif"))
##将结果导出
ggsave("nc.pdf",width = 8,height = 6)
ggsave("nc.tiff",width = 8,height = 6)
 ```
 ### 以上代码和例子数据在testdata目录
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
 
 
 
 
  
  
