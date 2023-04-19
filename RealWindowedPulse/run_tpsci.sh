#!/bin/bash

#SBATCH --nodes=1
#SBATCH -p normal_q
#SBATCH --cpus-per-task=16
#SBATCH -t 60:00:00
#SBATCH --account=nmayhall_group
##SBATCH --account=nmayhall_group-paid
#SBATCH --job-name p1_17
#SBATCH --exclusive

export NTHREAD=16
export JULIAENV=/home/nbraunsc/FermiCG/


# allow time to load modules, this takes time sometimes for some reason
sleep 10

hostname

module reset
module load Python/3.8.6-GCCcore-10.2.0
module load Julia/1.7.2-linux-x86_64

echo "Usage: sbatch submit.sh {input file} {data file} {data file}"

# set these to 1 if you are using julia threads heavily to avoid oversubscription
export MKL_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

export INFILE=$1
export OUTFILE="${INFILE}.out"
export WORKDIR=$(pwd)

echo $INFILE
echo $OUTFILE
echo $WORKDIR
echo $TMPDIR

cp $INFILE $TMPDIR/
if [ "$2" ]
then
        cp $2 $TMPDIR/
fi
if [ "$3" ]
then
        cp $3 $TMPDIR/
fi
cd $TMPDIR


#Start an rsync command which runs in the background and keeps a local version of the output file up to date
touch $OUTFILE
while true; do rsync -av $OUTFILE $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"; sleep 60; done &


julia --project=$JULIAENV -t $NTHREAD $INFILE >& $OUTFILE

cp $OUTFILE $WORKDIR/"${INFILE}.out"
rm $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.out"

mkdir $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

cp -r * $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"
rm -r *

#moving standard output slurm file to specific job directory
mv $WORKDIR/"slurm-${SLURM_JOB_ID}.out" $WORKDIR/"${INFILE}.${SLURM_JOB_ID}.scr"

exit

