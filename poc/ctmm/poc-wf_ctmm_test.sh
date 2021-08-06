#---- Test ctmm script locally ----#

wd=~/projects/ms3/analysis/full_workflow_poc/test4
src=~/projects/ms3/src

cd $wd

dat=data/seg_dates.csv
out='ctmm'

#Sequential
$src/poc/ctmm/poc_ctmm.r $dat $out -t

#Parallel
$src/poc/ctmm/poc_ctmm.r $dat $out -p mc -c 4 -t 

#--- Prep for hpc ----#

#Upload code to hpc
cd $src
git status
git add .
git commit -am 'add timing file output'
git push

#Upload project files
wdx="~/projects/ms3/analysis/full_workflow_poc/test4"

cd $wd

ssh grace "mkdir -p $wdx" #make sure to use double not single quotes!
scp -r ctfs grace:$wdx
scp -r data grace:$wdx

#--- HPC ---#

ssh grace
cd ~/projects/ms3/src
git pull

#--- Interactive queue

srun --pty -p interactive -n 4 bash #request four tasks in the interactive queue

wd=~/projects/ms3/analysis/full_workflow_poc/test4
src=~/projects/ms3/src

cd $wd

module load miniconda
source activate parallelR3

dat=data/seg_dates.csv
out='ctmm'
logs=logs

#Sequential execution
Rscript --vanilla $src/poc/ctmm/poc_ctmm.r $dat $out -t

#Parallel execution
mpirun -n 4 R --slave -f $src/poc/ctmm/poc_ctmm.r --args $dat $out -p mpi -m $logs -t

#--- Scavence queue
wd=~/projects/ms3/analysis/full_workflow_poc/test4

#need these for all scripts
export src=~/projects/ms3/src
export dat=data/seg_dates.csv
export out=ctmm
export logs=$out/mpilogs

#slurm variables
export n=4
export t=10:00
export p=scavenge
export J=test4
export mail=NONE

cd $wd

# These have to start with the --option b/c echo won't print - as first character
pars=`echo --ntasks $n -p $p -t $t -J $J --mail-type $mail`
exp=`echo --export=ALL,n=$n,p=$p,t=$t,J=$J,mail=$mail`

sbatch $pars $exp $src/poc/ctmm/poc_ctmm_sbatch.sh