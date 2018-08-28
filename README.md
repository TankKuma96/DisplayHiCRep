# DisplayHiCRep
R package for displaying result which got from r package 'Hicrep'

Tiankai Xie <xtk96kevin@gmail.com>

# Introduction

After tesing and establishing the feasibility of the method that were originally developed for measuring reproducibility of Hi-C data to HiChIP data, namely, HiCRep, which was developed by our lab(link: https://github.com/MonkeyLB/hicrep), we developed an expansion package, namely, `DisplayHiCRep`. It is used for calculating the reproducibility of HiChIP data and displaying all inter-chromosome and intra-chromosome in data frame and heat map. 

There are 2 functions wrapped inside, `juicer2mat` for HiChIP data with .hic data format and `cooler2mat` for HiChIP data with cool data format. Both functions will show the reproducibility scores for both intra-chromosome and inter-chromosome in data frame format. Plus, inside both functions, there is an opinion for controlling whether people need heat map based on results or not, to give out more choices.

# Citation

Cite the paper of our lab:

HiCRep: assessing the reproducibility of Hi-C data using a 
stratum-adjusted correlation coefficient. Tao Yang, Feipeng Zhang, Galip
Gurkan Yardimci, Fan Song, Ross C Hardison, William Stafford Noble, 
Feng Yue, Qunhua Li. Genome Research 2017. doi: 10.1101/gr.220640.117.

# Installation
Download the source package [DisplayHiCRep_0.0.0.9000.tar.gz](https://github.com/TankKuma96/DisplayHiCRep/blob/master/DisplayHiCRep_0.0.0.9000.tar.gz) from Github.

Or install from R with:

```{r}
library(devtools)
install_github("TankKuma96/DisplayHiCRep")
library(DisplayHiCRep)
```
 

# Data processing
1). Firstly, generating readable files by Juicer or Cooler on Linux based platform.

  -- for .hic files, use Juicer to dump all .txt files which contain all non-repeating  comparison data between inter-chromosome and intra-chromosome data.
  
  -- for .cool files, use Cooler to generate a huge .txt file which contain all non-repeating data about comparison between inter-chromosome and intra-chromosome.
  
2). Read files into right format

  -- for .hic files, use all other parameters provided to read all .txt files into same amount of data frames
  -- for .cool files, only take three columns of that huge files: chromosome 1 starting position, chromosome 2 starting position and counts. Then, group all the rows that contain same starting points of two chromosomes, and separate that huge file into  data frames that has same amount as the number of comparisons between inter-chromosome and intra-chromosome data.
  
3). Form huge super matrix by using all those data frames just created.

4). Use `fast.mean.filter` in HiCRep to smooth that super matrix.

5). Split the smoothed super matrix into data frames again.

6). Calculate reproducibility scores by 2 groups : diagonal group and off-diagonal group. 

  -- for diagonal group, use `get.scc` in HiCRep to calculate its SCC
  
  -- for off-diagonal group, use `cor` and pick option: Pearson correlation to get its correlations.

Figure1. `DisplayHiCRep` data processing
                          
![Figure1. `DisplayHiCRep` data processing](https://github.com/TankKuma96/DisplayHiCRep/blob/master/vignettes/new_dp.png)
