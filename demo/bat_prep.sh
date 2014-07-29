home=..
data=$home/data
demo=$home/demo
#domain=aimpeak
domain=hadsst
#domain=sarcos
seed=5
mode=1
blk=2
blksize=1000
sptset=32
percent=10
./prep -dataset $data/$domain -output $demo/$domain -seed $seed -mode $mode -machine $blk -blocksize $blksize -percent $percent -support $sptset
