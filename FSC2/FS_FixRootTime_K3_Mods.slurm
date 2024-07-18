#!/bin/bash
#SBATCH --job-name=fastsimcoal.noreponly_v2
#SBATCH --partition=kill-shared
#SBATCH --time=03-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=NONE
#SBATCH --mail-user=hwsung@hawaii.edu
#SBATCH --error=job%A_%a.err ## %A - filled with jobid
#SBATCH --output=job%A_%a.out ## %A - filled with jobid
#SBATCH --array=1-1250

source ~/.bash_profile

cd /home/hwsung/kantar_koastore/helen/analyses/fastsimcoal/noreponly_v2/Reps  # set the correct working directory - make sure that the array number is correct, too


for x in */*; do
  repdirs=(${repdirs[@]} "${x}")
done

cd ${repdirs[($SLURM_ARRAY_TASK_ID-1)]}   ## here, I subtract 1 from the $SLURM_ARRAY_TASK_ID because bash indexing starts at zero - I think it's less confusing to subtract 1 here than to remember to do it when specifying the number of jobs for the array

est=$(echo *.est)
tpl=$(echo *.tpl)

fsc26 -t $tpl -e $est -n 100000 -L 40 -m -M -u -x  -q 