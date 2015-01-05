#!/bin/sh
######################################################################################
##                  This script demos how to use the PGPR code.
##   To run the PGPR code, there are some common setting for the experiments, including 
##   the working path, parameters for the algorithms, and the configuration of MPICH,
## 
##   For the format specification of the original input, please refer to the readme.txt in
##   the ../data folder. Before calling the main approximated gp algorithms, we need to 
##   call ./prep to generate the train dataset (.trn), test dataset (.tst),  support dataset
##   (.spt if mode=1), and hyperparameter file (.hyp).
## 
##   Each GP algorithm would take generated files by ./prep as the input file. Except for 
##   displaying the general output information (incuring time, root mean square error (RMSE),
##   mean negative log probability (MNLP)), the gp algorithm will also generate a .rts file, 
##   each line of which has two columns, mean and variance for the corresponding test points.
##   

######### Basic setting for the experiments ################
###  Below are the options for the experiments setting, you may need to change them for your own
###  experiments setting. 
###      data: the path for the data file
###      demo: the path for the demo folder
###      mode: since the support set selecting process would be slow if support set is large,
###           to avoid regenerating the support set, you can set mode = 0
###       dom: the domain name of the using data.
###       blk: the number of blocks
###   blksize: the size of each block
### bandwidth: the bandwidth for LMA
###   percent: the percent of total data used for the testing
###      rank: the approximated rank 
###      sset: the number of support data points
###     hosts: the configure of the hosts for the MPICH

home=..
data=$home/data                
demo=$home/demo               
#dom=aimpeak
#dom=hadsst
dom=sarcos                    
seed=5
mode=1                        
blk=2                         
blksize=100
bandwidth=0
percent=10                    
hyp=$dom.hyp
inf=$dom.5
rank=512                       
sset=32    
hosts=hosts 
host=$(head -n1 $hosts)

#prepare the train data, test set and support set files
echo ===Prepare the train data, test set and support set files===
echo Domain: $dom
echo Machines/Blocks: $blk
echo Blocksize: $blksize
echo Bandwidth: $bandwidth
echo Size of support set: $sset
echo Reduced rank: $rank
echo Percentage of data for test: $percent
echo Random seed: $seed
./prep -dataset $data/$dom -output $demo/$dom -seed $seed -mode $mode -machine $blk -blocksize $blksize -percent $percent -support $sset

echo ===============Run Gaussian process regressions=============
#Run parallel Gaussian process regressions and their sequential versions
echo " ALGO | Runtime |  RMSE  |  MNLP  |"
./pic -hyper $hyp -in $inf -out pic_$inf  -blocks $blk
mpirun -f $hosts -n $blk ./ppic -hyper $hyp -in $inf -out ppic_$inf -blocks $blk
./pitc -hyper $hyp -in $inf -out pitc_$inf -blocks $blk
mpirun -f $hosts -n $blk ./ppitc -hyper $hyp -in $inf -out ppitc_$inf -blocks $blk
mpirun -f $hosts -n $blk ./plma -hyper $hyp -in $inf -out plma_$inf -blocks $blk -bandwidth $bandwidth
#Run full Gaussian process regression. If the training data is to large
#(>2000), it's better to comment out the FGP which requires extremely
#long time.  
./fgp -hyper $hyp -in $inf -out fgp_$inf 

